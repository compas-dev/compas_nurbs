# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
import sphinx_compas_theme

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]
if os.getenv('SPELLCHECK'):
    extensions += 'sphinxcontrib.spelling',
    spelling_show_suggestions = True
    spelling_lang = 'en_US'

source_suffix = '.rst'
master_doc = 'index'
project = 'compas_nurbs'
year = '2018'
author = 'Gramazio Kohler Research'
copyright = '{0}, {1}'.format(year, author)
version = release = '0.1.0'

pygments_style = 'trac'  # Perhaps change to sphinx
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/gramaziokohler/compas_nurbs/issues/%s', '#'),
    'pr': ('https://github.com/gramaziokohler/compas_nurbs/pull/%s', 'PR #'),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
html_theme = 'compas'
html_theme_path = sphinx_compas_theme.get_html_theme_path()
html_theme_options = {
    'logo': 'logo.png',
    'description': 'compas_nurbs',
    'github_user': 'gramaziokohler',
    'github_repo': project,
    'fixed_sidebar': True,
}

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_static_path = ['_static']
html_sidebars = {
    '**': ['about.html', 'navigation.html', 'searchbox.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
