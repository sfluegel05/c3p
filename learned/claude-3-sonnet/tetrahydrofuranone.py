"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:33524 tetrahydrofuranone
Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for tetrahydrofuran ring with an oxo substituent
    oxo_thf_pattern = Chem.MolFromSmarts("[OX1]1[C@@H][C@@H][C@@H][C@@H]O1")
    thf_ring_atoms = mol.GetSubstructMatches(oxo_thf_pattern)
    
    # Check if there are any additional rings
    additional_rings = Chem.GetSymmSSSR(mol)
    additional_rings = [ring for ring in additional_rings if set(ring) != set(thf_ring_atoms[0])]

    if not thf_ring_atoms:
        return False, "No tetrahydrofuran ring with oxo substituent found"
    elif len(additional_rings) > 1:
        return False, "More than one additional ring system present"

    # Check for specific substituent patterns
    substituent_pattern = Chem.MolFromSmarts("[OX2]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    
    # Filter out substituents that are part of the tetrahydrofuran ring
    substituent_matches = [match for match in substituent_matches if set(match) != set(thf_ring_atoms[0])]

    # Check for common functional groups
    ethers = Chem.MolFromSmarts("[OD2]([#6])[#6]")
    esters = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    alcohols = Chem.MolFromSmarts("[OX2H]")
    
    has_ether = mol.HasSubstructMatch(ethers)
    has_ester = mol.HasSubstructMatch(esters)
    has_alcohol = mol.HasSubstructMatch(alcohols)

    # Check ring size and heteroatom count
    thf_ring = mol.GetAtomRingInfo().IsCyclic()
    heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)

    if not thf_ring or len(thf_ring_atoms[0]) != 5 or heteroatom_count > 2:
        return False, "Ring size or heteroatom count does not match tetrahydrofuranone"

    return True, "Molecule contains a tetrahydrofuran ring with an oxo substituent"