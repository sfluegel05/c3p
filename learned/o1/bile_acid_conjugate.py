"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a hydrophilic group such as glycine,
    taurine, amino acids, sulfate, glucuronic acid, glucose, other sugars, or coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
        
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a general pattern for the steroid nucleus (four fused rings)
    steroid_pattern_smarts = '[#6]12CC[C@H]3[C@@H](C1)CCC4=C3C=CC=C4C2'
    steroid_pattern = Chem.MolFromSmarts(steroid_pattern_smarts)
    if steroid_pattern is None:
        return False, "Invalid steroid nucleus SMARTS pattern"
    
    # Check for steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"
    
    # Define a pattern for the bile acid side chain ending with carboxylic acid or derivative
    side_chain_smarts = '[C@H](CCC(=O)[O;H,-])[C@@H](C)CC(=O)[O;H,-]'
    side_chain_pattern = Chem.MolFromSmarts(side_chain_smarts)
    if side_chain_pattern is None:
        return False, "Invalid side chain SMARTS pattern"
    
    # Combine steroid nucleus and side chain to define bile acid core
    bile_acid_smarts = steroid_pattern_smarts + side_chain_smarts
    bile_acid_pattern = Chem.MolFromSmarts(bile_acid_smarts)
    if bile_acid_pattern is None:
        return False, "Invalid bile acid core SMARTS pattern"
    
    # Check for bile acid core
    if not mol.HasSubstructMatch(bile_acid_pattern):
        return False, "No bile acid core found"

    # Define SMARTS patterns for conjugated groups attached via amide or ester linkage
    conjugates = {
        "glycine": 'C(=O)NCC(=O)[O;H,-]',
        "taurine": 'C(=O)NCCS(=O)(=O)[O;H,-]',
        "amino_acid": 'C(=O)N[C@@H]([#6])C(=O)[O;H,-]',  # General amino acid
        "sulfate": 'OS(=O)(=O)[O;H,-]',
        "glucuronic_acid": '[#6]-1(-[#8]-[#6]-[#8]C(=O)[O;H,-])-[#6]([#8])-[#6]([#8])-[#6]([#8])-[#6]-1',
        "glucose": '[#6]-1(-[#8]-[#6])-[#6]([#8])-[#6]([#8])-[#6]([#8])-[#6]-1',
        "coenzyme_A": 'NC(=O)CCNC(=O)CC(=O)NCCS',
    }

    # Check for conjugation with any of the specified groups
    conjugated = False
    conjugate_name = ""
    for name, smarts in conjugates.items():
        conjugate_pattern = Chem.MolFromSmarts(smarts)
        if conjugate_pattern is None:
            continue  # Skip invalid pattern
        # Attachment via an amide or ester bond to the bile acid carboxyl group
        linker_smarts = 'C(=O)[O,N]-' + smarts
        linker_pattern = Chem.MolFromSmarts(linker_smarts)
        if linker_pattern is None:
            continue  # Skip invalid linker pattern
        if mol.HasSubstructMatch(linker_pattern):
            conjugated = True
            conjugate_name = name
            break

    if not conjugated:
        return False, "No conjugated hydrophilic group found"

    return True, f"Contains bile acid core conjugated with {conjugate_name}"