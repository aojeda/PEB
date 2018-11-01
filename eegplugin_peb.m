% eegplugin_peb() - Plugin for solving the inverse problem using the extended Parametric Empirical Bayes algorithm (PEB+).
% Usage:
%   >> eegplugin_peb(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Alejandro Ojeda, NEATLabs, UCSD, 2018
%
% See also: eeglab()

function vers = eegplugin_peb(fig,try_strings, catch_strings)
vers = 'PEB+1.0.1';
p = fileparts(which('PEB'));
addpath(genpath(p));
h = findobj(gcf, 'tag', 'tools');
hmMenu = uimenu( h, 'label', 'PEB+');
uimenu( hmMenu, 'label', 'EEG source estimation','callback','EEG = pop_inverseSolution(EEG);');
