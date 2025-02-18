"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:73404 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonic acid group (-S(=O)(=O)O)
    sulfonic_acid = Chem.MolFromSmarts("[S](=[O])(=[O])[O]")
    if not mol.HasSubstructMatch(sulfonic_acid):
        return False, "No sulfonic acid group found"

    # Verify sulfonic acid is attached to carbon via C-S bond
    sulfur_atoms = [atom.GetIdx() for atom in mol.GetAtoms() 
                   if atom.GetAtomicNum() == 16 and 
                   any(bond.GetBondType() == Chem.BondType.SINGLE 
                       for bond in atom.GetBonds())]
    
    valid_sulfur = False
    for s_idx in sulfur_atoms:
        sulfur = mol.GetAtomWithIdx(s_idx)
        # Check if sulfur is bonded to carbon
        if any(neighbor.GetAtomicNum() == 6 for neighbor in sulfur.GetNeighbors()):
            valid_sulfur = True
            break
    
    if not valid_sulfur:
        return False, "Sulfonic acid not attached to carbon via C-S bond"

    # Check lipid characteristics (long hydrocarbon chains)
    # Count total carbons and look for long chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Insufficient carbons ({c_count}) for lipid structure"

    # Check for lipid-like features (long aliphatic chains)
    # Look for at least one chain of 8+ carbons
    long_chain = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long hydrocarbon chain found (minimum 8 carbons)"

    # Check for lipid-like functional groups (esters/amides)
    ester = Chem.MolFromSmarts("[O][C]=O")
    amide = Chem.MolFromSmarts("[N][C]=O")
    if not mol.HasSubstructMatch(ester) and not mol.HasSubstructMatch(amide):
        return False, "No ester/amide groups characteristic of lipids"

    return True, "Contains sulfonic acid group attached via C-S bond to lipid structure"