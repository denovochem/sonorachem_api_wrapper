import random
import re
from rdkit import Chem
from rdkit import RDLogger
import logging
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog('rdApp.*')
logging.disable(logging.CRITICAL)

def MolWithoutIsotopesToSmiles(mol):
  atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]
  for atom, isotope in atom_data:
      if isotope:
          atom.SetIsotope(0)
  return mol


def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    clean_mol = rdMolStandardize.Cleanup(mol) 
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
    mol =  MolWithoutIsotopesToSmiles(uncharged_parent_clean_mol)
    return mol

def randomize_smiles(smiles, standardize_mol=False, isomeric=True, shuffle_order=True, remove_mapping=True):
    try:
        x = smiles.split('.')
        if shuffle_order:
            random.shuffle(x)
        frags = []
        for i in x:
            if standardize_mol:
                m = standardize(i)
            else:
                m = Chem.MolFromSmiles(i)
            if remove_mapping:
                [a.SetAtomMapNum(0) for a in m.GetAtoms()]
            new_atom_order = list(range(m.GetNumAtoms()))
            random.shuffle(new_atom_order)
            random_mol = Chem.RenumberAtoms(m, newOrder=new_atom_order)
            random_smiles_string = str(Chem.MolToSmiles(random_mol, canonical=False, isomericSmiles=isomeric))
            frags.append(random_smiles_string)
        random_smiles_string = '.'.join(i for i in frags)
        return(random_smiles_string)
    except:
        pass

def canonicalize_smiles(smiles, standardize_mol=False, isomeric=True, remove_mapping=True):
    try:
        x = smiles.split('.')
        x = sorted(x)
        frags = []
        for i in x:
            if standardize_mol:
                m = standardize(i)
            else:
                m = Chem.MolFromSmiles(i)
            if remove_mapping:
                [a.SetAtomMapNum(0) for a in m.GetAtoms()]
            canonical_smiles_string = str(Chem.MolToSmiles(m, canonical=True, isomericSmiles=isomeric))
            frags.append(canonical_smiles_string)
        canonical_smiles_string = '.'.join(i for i in sorted(frags))
        return(canonical_smiles_string)
    except:
        pass

def randomize_reaction_smiles(smiles, standardize_mol=False, isomeric=True, grouped=False, shuffle_order=True):
    try:
        split_roles = smiles.split('>>')
        if grouped:
            randomized_rxn = '{'
        else:
            randomized_rxn = ''
        for x in split_roles:
            if x != '':
                y = x.split('}.{')
                if shuffle_order:
                    random.shuffle(y)
                for z in y:
                    if grouped:
                        z = z.replace('}', '').replace('{', '')
                    randomized_smiles = randomize_smiles(z, standardize_mol, isomeric)
                    randomized_rxn += randomized_smiles
                    if grouped:
                        randomized_rxn += '}.{'
                    else:
                        randomized_rxn += '.'
                if grouped:
                    randomized_rxn = randomized_rxn[:-2]
                    randomized_rxn += '>>{'
                else:
                    randomized_rxn = randomized_rxn[:-1]
                    randomized_rxn += '>>'
            else:
                randomized_rxn = randomized_rxn.rstrip('{')
                randomized_rxn += '>>'
        if grouped:
            randomized_rxn = randomized_rxn[:-3]
        else:
            randomized_rxn = randomized_rxn[:-2]
        return(randomized_rxn)
    except:
        pass

def canonicalize_reaction_smiles(smiles, standardize_mol=False, isomeric=True, grouped=False):
    try:
        split_roles = smiles.split('>>')
        if grouped:
            canonical_rxn = '{'
        else:
            canonical_rxn = ''
        for x in split_roles:
            if x != '':
                y = x.split('}.{')
                y = sorted(y)
                for z in y:
                    if grouped:
                        z = z.replace('}', '').replace('{', '')
                    canonical_smiles = canonicalize_smiles(z, standardize_mol, isomeric)
                    canonical_rxn += canonical_smiles
                    if grouped:
                        canonical_rxn += '}.{'
                    else:
                        canonical_rxn += '.'
                if grouped:
                    canonical_rxn = canonical_rxn[:-2]
                    canonical_rxn += '>>{'
                else:
                    canonical_rxn = canonical_rxn[:-1]
                    canonical_rxn += '>>'
            else:
                canonical_rxn = canonical_rxn.rstrip('{')
                canonical_rxn += '>>'
        if grouped:
            canonical_rxn = canonical_rxn[:-3]
        else:
            canonical_rxn = canonical_rxn[:-2]
        return(canonical_rxn)
    except:
        pass
