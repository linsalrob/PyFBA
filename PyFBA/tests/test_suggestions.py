import os
import unittest
import PyFBA

test_file_loc = os.path.dirname(os.path.realpath(__file__))


class SuggestionTest(unittest.TestCase):
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data('gramnegative', verbose=False)
    compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(organism_type='gramnegative')

    def test_essential(self):
        """Test the essential reactions"""
        suggested = PyFBA.gapfill.suggest_essential_reactions()
        self.assertEqual(len(suggested), 109)

    def test_limit_by_compound(self):
        """Test limiting adding reactions by the compounds that are present"""
        suggested = PyFBA.gapfill.suggest_by_compound(self.__class__.modeldata,
                                                      {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22712)
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
        toy_media = {PyFBA.metabolism.Compound('cpd00013', 'NH3', 'e')}
        suggested = PyFBA.gapfill.suggest_from_media(self.__class__.modeldata, {}, toy_media)
        self.assertEqual(len(suggested), 19)

    def test_orphan_compound(self):
        """Test suggestions based on orphan compounds"""
        suggested = PyFBA.gapfill.suggest_by_compound(self.__class__.modeldata,
                                                      {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22712)

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
        self.assertGreaterEqual(len(suggested), 2766)

    def test_roles_from_file(self):
        """Test suggestions based on roles in a file"""
        self.assertTrue(os.path.exists(os.path.join(test_file_loc, 'roles.txt')))
        suggs = PyFBA.gapfill.suggest_from_roles(os.path.join(test_file_loc, 'roles.txt'), self.__class__.reactions)
        self.assertEqual(len(suggs), 3)
        for rxn in {'rxn00577', 'rxn01192', 'rxn00799'}:
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
