# ----------------------------------------------------------------------------
# Copyright (c) 2026, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase


class ExamplesTests(TestPluginBase):
    """Run every usage example registered on the plugin's actions."""
    package = 'q2_confluence.tests'

    def test_examples(self):
        self.execute_examples()
