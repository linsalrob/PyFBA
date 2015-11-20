import os
import unittest

from gapfill import orphan_compound, essentials, limit_reactions, probability, suggest_reactions_without_proteins, \
    suggest_reactions_with_proteins, media, roles, subsystem
from metabolism import Compound
from parse import model_seed


class SuggestionTest(unittest.TestCase):

    compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()

    def test_essential(self):
        suggested = essentials.suggest_essential_reactions()
        self.assertEqual(len(suggested), 109)

    def test_limit_by_compound(self):
        suggested = orphan_compound.suggest_by_compound(self.__class__.compounds, self.__class__.reactions,
                                                        {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)
        limited = limit_reactions.limit_reactions_by_compound(self.__class__.reactions, {'rxn00001'}, suggested, 5)
        self.assertGreaterEqual(len(limited), 22712)

    def test_with_proteins(self):
        suggested = suggest_reactions_with_proteins(self.__class__.reactions)
        self.assertGreaterEqual(len(suggested), 2539)

    def test_without_proteins(self):
        suggested = suggest_reactions_without_proteins(self.__class__.reactions)
        self.assertGreaterEqual(len(suggested), 32157)

    def test_media(self):
        toy_media = {Compound('NH3', 'e')}
        suggested = media.suggest_from_media(self.__class__.compounds, {}, toy_media)
        self.assertEqual(len(suggested), 758)

    def test_orphan_compound(self):
        suggested = orphan_compound.suggest_by_compound(self.__class__.compounds, self.__class__.reactions,
                                                        {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)

    def test_probability(self):
        self.assertTrue(os.path.exists('tests/reaction_list.txt'))
        reactions2run = set()
        with open('tests/reaction_list.txt', 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                if "biomass" in l.lower():
                    continue
                r = l.strip()
                if r in self.__class__.reactions:
                    reactions2run.add(r)
        suggested = probability.compound_probability(self.__class__.reactions, reactions2run)
        self.assertEqual(len(suggested), 727)

    def test_roles_from_file(self):
        self.assertTrue(os.path.exists('tests/roles.txt'))
        suggestions = roles.suggest_from_roles('tests/roles.txt', self.__class__.reactions)
        self.assertEqual(len(suggestions), 4)
        for rxn in {'rxn00799', 'rxn01192', 'rxn00577', 'rxn01277'}:
            self.assertIn(rxn, suggestions)

        suggestions = roles.suggest_from_roles('tests/roles.txt', self.__class__.reactions, 0.9)
        self.assertEqual(len(suggestions), 1)
        for rxn in {'rxn01277'}:
            self.assertIn(rxn, suggestions)

    def test_subsystems(self):
        r2r = {'rxn00006', 'rxn00189', 'rxn00194', 'rxn00206', 'rxn00322', 'rxn00405', 'rxn05216', 'rxn10136',
               'rxn17196', 'rxn20595', 'rxn26755'}
        suggestions = subsystem.suggest_reactions_from_subsystems(self.__class__.reactions, r2r)
        self.assertGreaterEqual(len(suggestions), 407)
        # not testing all 407 reactions!
        for r in ['rxn00006', 'rxn00206', 'rxn17196']:
            self.assertIn(r, suggestions)
        suggestions = subsystem.suggest_reactions_from_subsystems(self.__class__.reactions, r2r, threshold=1)
        self.assertEqual(len(suggestions), 3)
        for r in ['rxn17196', 'rxn26755', 'rxn20595']:
            self.assertIn(r, suggestions)

if __name__ == '__main__':
    unittest.main()
