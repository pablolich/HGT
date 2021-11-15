#!/usr/bin/env python3

__appname__ = '[assembly_cr.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys

## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #1. Initialize random community. Let it evolve towards a globally stable 
    #equilibrium by prunning the system through dynamics to a feasible state. 
    #(for this system feasiblity ==> global stability)
    #Size of the community 
    m = 10
    #2. Mutate invader, check for invasibility. 
    #3. (provided it can invade) Let the system evolve towards a new global
    #attractor
    #4. spread mutation throughout community by changing C with matrix H
    #5. Let the system evolve towards another equilibrium
    #6. Repeat until convergence of some sort is reached. Maybe this
    #this convergence is when after n_max invasion attempts, there are no 
    #successes. The community would have reached an end point. 

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

