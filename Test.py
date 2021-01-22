# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 17:02:25 2021

@author: Silvia Vargas
"""
from Diffusion_2D import Diffusion
import sys
from sys import argv
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

dx = config.get('parameters', 'dx')
nu = config.get('parameters', 'nu')
kind = config.get('parameters', 'kind')
nt = config.get ('parameters', 'nt')
L = config.get ('parameters', 'L')

dx = float(dx)
nu = float(nu)
nt = int(nt)
L = float(L)

Diff1=Diffusion(dx,nu,kind,nt,L)
Initial_conditions = Diff1.u[:,:,0]
Diff1.evolve_ts()
Final_timestep = Diff1.u[:,:,Diff1.nt]

def test_Initial_conditions_are_square(Initial_conditions):
    """ 
    Check that the Initial conditions are a square matrix, and therefore dx and dy are equal
    """
    assert all(len(row) == len(Initial_conditions) for row in Initial_conditions)


def test_Initial_conditions_contain_only_zeros_and_ones(Initial_conditions):
    """
    Check that the Initial conditions of the class Diffusion contain only 0s and 1s
    """
     # flatten up the mesh into one single list
     # and set on the list it should be [0,1] if it
     # contains only 0 and 1. Then do sum on that will 
     # return 1
     
    assert (sum(set(sum(Initial_conditions,[])))  <= 1)
    
def test_Initial_conditions_are_the_first_to_appear(Initial_conditions):
    """Check that the initial condition has index 0"""
    
def test_amount_of_matter_is_the_same(Initial_conditions, Final_timestep):
    """
    Test if the amount of matter is the same before and after diffusion
    """
    assert round((sum(sum(Initial_conditions)))) == round((sum(sum(Final_timestep))))
    
