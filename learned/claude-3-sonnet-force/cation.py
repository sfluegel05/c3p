"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:49808 cation
A monoatomic or polyatomic species having one or more elementary charges of the proton.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check formal charge
    formal_charge = AllChem.GetFormalCharge(mol)
    if formal_charge > 0:
        return True, f"Molecule has positive formal charge of {formal_charge}"

    # Check for common cationic functional groups
    cationic_patterns = (
        "[N+]",       # Quaternary ammonium ions
        "[N+]#N",     # Diazonium ions
        "[NH3+]",     # Ammonium ions
        "[O+]",       # Oxonium ions
        "[S+]",       # Sulfonium ions
        "[P+]",       # Phosphonium ions
        "[Se+]",      # Selenonium ions
        "[Te+]",      # Telluronium ions
        "[As+]",      # Arsonium ions
        "[Cl+]",      # Chloronium ions
        "[I+]",       # Iodonium ions
        "[Br+]",      # Bromonium ions
    )
    for pattern in cationic_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return True, f"Molecule contains {pattern} cationic group"

    # Check for zwitterionic species
    zwitterion_pattern = Chem.MolFromSmarts("[+,-]")
    if mol.HasSubstructMatch(zwitterion_pattern):
        return True, "Molecule is a zwitterionic species"

    # Check for aromatic heterocyclic cations
    aromatic_cations_pattern = Chem.MolFromSmarts("[n+,o+,s+,se+,te+]")
    if mol.HasSubstructMatch(aromatic_cations_pattern):
        return True, "Molecule contains aromatic heterocyclic cation"

    # Check for resonance-stabilized cations
    resonance_cations_pattern = Chem.MolFromSmarts("[c,n,o,s,se,te][c+]")
    if mol.HasSubstructMatch(resonance_cations_pattern):
        return True, "Molecule contains resonance-stabilized cation"

    # If no cationic features found, classify as non-cation
    return False, "No cationic features found"