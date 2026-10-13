# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_confluence.plugin_setup import plugin


class SignatureTests(unittest.TestCase):

    def test_output_names_distinct_from_inputs_and_parameters(self):
        # TODO: this is AI-generated code; human review pending
        # q2cli classifies each signature entry by name, so an output named
        # like an input or parameter is rendered as that input and gets no
        # --o- option (frequency's output was named table, like its input).
        for name, action in plugin.actions.items():
            sig = action.signature
            taken = set(sig.inputs) | set(sig.parameters)
            clash = taken & set(sig.outputs)
            self.assertEqual(clash, set(), f'{name}: {clash}')


if __name__ == '__main__':
    unittest.main()
