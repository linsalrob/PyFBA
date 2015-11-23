import os
import unittest
import PyFBA

test_file_loc = ''
if os.path.exists('tests/roles.txt'):
    test_file_loc = 'tests'
elif os.path.exists('PyFBA/tests/roles.txt'):
    test_file_loc = 'PyFBA/tests'


class SuggestionTest(unittest.TestCase):
    compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes()

    def test_essential(self):
        """Test the essential reactions"""
        suggested = PyFBA.gapfill.suggest_essential_reactions()
        self.assertEqual(len(suggested), 109)

    def test_limit_by_compound(self):
        """Test limiting adding reactions by the compounds that are present"""
        suggested = PyFBA.gapfill.suggest_by_compound(self.__class__.compounds, self.__class__.reactions,
                                                      {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)
        limited = PyFBA.gapfill.limit_reactions_by_compound(self.__class__.reactions, {'rxn00001'}, suggested, 5)
        self.assertGreaterEqual(len(limited), 22712)

    def test_with_proteins(self):
        """Test suggesting reactions with proteins"""
        suggested = PyFBA.gapfill.suggest_reactions_with_proteins(self.__class__.reactions)
        self.assertGreaterEqual(len(suggested), 2539)

    def test_without_proteins(self):
        """Test suggesting reactions without proteins"""
        suggested = PyFBA.gapfill.suggest_reactions_without_proteins(self.__class__.reactions)
        self.assertGreaterEqual(len(suggested), 32157)

    def test_media(self):
        """Test suggestions based on a media"""
        toy_media = {PyFBA.metabolism.Compound('NH3', 'e')}
        suggested = PyFBA.gapfill.suggest_from_media(self.__class__.compounds, self.__class__.reactions, {}, toy_media)
        self.assertEqual(len(suggested), 758)

    def test_orphan_compound(self):
        """Test suggestions based on orphan compounds"""
        suggested = PyFBA.gapfill.suggest_by_compound(self.__class__.compounds, self.__class__.reactions,
                                                      {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)

    def test_probability(self):
        """Test suggestions based on probability of the reactions occurring"""
        self.assertTrue(os.path.exists(os.path.join(test_file_loc, 'reaction_list.txt')))
        reactions2run = set()
        with open(os.path.join(test_file_loc, 'reaction_list.txt'), 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                if "biomass" in l.lower():
                    continue
                r = l.strip()
                if r in self.__class__.reactions:
                    reactions2run.add(r)
        suggested = PyFBA.gapfill.compound_probability(self.__class__.reactions, reactions2run)
        self.assertEqual(len(suggested), 727)

    def test_roles_from_file(self):
        """Test suggestions based on roles in a file"""
        self.assertTrue(os.path.exists(os.path.join(test_file_loc, 'roles.txt')))
        suggs = PyFBA.gapfill.suggest_from_roles(os.path.join(test_file_loc, 'roles.txt'), self.__class__.reactions)
        self.assertEqual(len(suggs), 4)
        for rxn in {'rxn00799', 'rxn01192', 'rxn00577', 'rxn01277'}:
            self.assertIn(rxn, suggs)

        suggs = PyFBA.gapfill.suggest_from_roles(os.path.join(test_file_loc, 'roles.txt'),
                                                 self.__class__.reactions, 0.9)
        self.assertEqual(len(suggs), 1)
        for rxn in {'rxn01277'}:
            self.assertIn(rxn, suggs)

    def test_subsystems(self):
        """Test suggestions based on subsystems that are covered"""
        r2r = {'rxn00006', 'rxn00189', 'rxn00194', 'rxn00206', 'rxn00322', 'rxn00405', 'rxn05216', 'rxn10136',
               'rxn17196', 'rxn20595', 'rxn26755'}
        suggestions = PyFBA.gapfill.suggest_reactions_from_subsystems(self.__class__.reactions, r2r)
        self.assertGreaterEqual(len(suggestions), 396)
        suggestions = PyFBA.gapfill.suggest_reactions_from_subsystems(self.__class__.reactions, r2r, threshold=0.3)
        self.assertEqual(len(suggestions), 2)


if __name__ == '__main__':
    unittest.main()
