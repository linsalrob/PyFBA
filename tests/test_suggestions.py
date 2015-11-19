import os
import unittest

from gapfill import orphan_compound, essentials, limit_reactions, probability, suggest_reactions_without_proteins, \
    suggest_reactions_with_proteins, media
from metabolism import Compound
from parse import model_seed


class SuggestionTest(unittest.TestCase):

    def test_orphan_compound(self):
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        suggested = orphan_compound.suggest_by_compound(compounds, reactions, {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)

    def test_essential(self):
        suggested = essentials.suggest_essential_reactions()
        self.assertEqual(len(suggested), 109)

    def test_limit_by_compound(self):
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        suggested = orphan_compound.suggest_by_compound(compounds, reactions, {'rxn00001'}, 2)
        self.assertGreaterEqual(len(suggested), 22729)
        limited = limit_reactions.limit_reactions_by_compound(reactions, {'rxn00001'}, suggested, 5)
        self.assertGreaterEqual(len(limited), 22712)

    def test_with_proteins(self):
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        suggested = suggest_reactions_with_proteins(reactions)
        self.assertGreaterEqual(len(suggested), 2539)

    def test_without_proteins(self):
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        suggested = suggest_reactions_without_proteins(reactions)
        self.assertGreaterEqual(len(suggested), 32157)

    def test_probability(self):
        self.assertTrue(os.path.exists('tests/reaction_list.txt'))
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        reactions2run = set()
        with open('tests/reaction_list.txt', 'r') as f:
            for l in f:
                if l.startswith('#'):
                    continue
                if "biomass" in l.lower():
                    continue
                r = l.strip()
                if r in reactions:
                    reactions2run.add(r)
        suggested = probability.compound_probability(reactions, reactions2run)
        self.assertEqual(len(suggested), 727)

    def test_media(self):
        compounds, reactions, enzymes = model_seed.compounds_reactions_enzymes()
        toymedia = {Compound('NH3', 'e')}
        suggested = media.suggest_from_media(compounds, {}, toymedia)
        self.assertEqual(len(suggested), 758)

if __name__ == '__main__':
    unittest.main()
