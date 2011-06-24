# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Clear memory and close all open plot windows
clear all;
close all;

# Ensure that the Sundials toolbox is in the path
if (~ exist ('CVode'))
  if (exist ('startup_STB'))
    startup_STB;
  else
    error ('%s: requires the Sundials toolbox, please see the README file for installation instructions');
  end
end