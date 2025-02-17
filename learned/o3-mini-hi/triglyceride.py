"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: Triglycerides
Definition: Any glyceride resulting from the condensation of all three hydroxy groups of glycerol (propane-1,2,3-triol) with fatty acids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined as a glyceride where all three hydroxy groups of glycerol are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a triglyceride, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Check for the presence of three ester groups.
    # Ester groups are represented by an oxygen single-bonded to a carbonyl group.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups; expecting exactly 3 for a triglyceride."
    
    # Step 2: Look for a glycerol backbone.
    # A glycerol backbone (propane-1,2,3-triol) in its esterified form usually shows a three-carbon chain
    # where each terminal carbon (primary carbons) and the middle carbon (secondary carbon) are connected to an ester group.
    # We use a simplified SMARTS pattern for a three-carbon chain.
    glycerol_pattern = Chem.MolFromSmarts("[$([CH2]-[CH]-[CH2])]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No suitable glycerol backbone (three-carbon chain) found."
    
    # Step 3: Check for the presence of fatty acid chains.
    # We require that each ester group leads to a relatively long aliphatic chain.
    # Here we use a simple pattern to find a chain of at least four carbons:
    fatty_chain_pattern = Chem.MolFromSmarts("[CX4]{4,}")
    fatty_chain_matches = mol.GetSubstructMatches(fatty_chain_pattern)
    if len(fatty_chain_matches) < 3:
        return False, f"Not enough fatty acid chains detected; found {len(fatty_chain_matches)} chain(s) with at least four carbons."
    
    # Step 4: Check molecule properties.
    # Typically, triglycerides have a fairly high number of rotatable bonds and a molecular weight above 500 Da.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds, suggesting fatty acid chains may be too short."
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a triglyceride."
    
    # Step 5: (Optional) Check for approximate counts of carbon and oxygen atoms.
    # A typical triglyceride should have three ester carbonyls (each with one oxygen) plus the oxygens in the glycerol backbone.
    # So altogether one expects roughly 6 oxygens.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for a triglyceride."
    
    # Count carbons (a TG should have long fatty acid chains, so expect many carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms for a triglyceride."

    # Passed all checks
    return True, "Molecule is classified as a triglyceride: it has a glycerol backbone with 3 esterified fatty acid chains."
    
# Example usage:
# print(is_triglyceride("O([C@H](COC(=O)CCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCCCCCCC"))