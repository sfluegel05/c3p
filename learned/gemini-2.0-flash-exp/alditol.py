"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all atoms
    atoms = mol.GetAtoms()

    #Find the linear backbone:
    chain_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CHX4](O)[CH2X4]") # minimum 4 carbons
    chain_matches = mol.GetSubstructMatches(chain_pattern)

    if not chain_matches:
        return False, "Molecule does not have a linear chain of carbons with OH groups"
    
    # We assume that the first match contains the alditol backbone
    chain_atoms = chain_matches[0]
    backbone_atoms = []
    for atom_idx in chain_atoms:
        backbone_atoms.append(mol.GetAtomWithIdx(atom_idx))
    
    #Check if the chain has CH2OH termini
    if not (backbone_atoms[0].GetTotalDegree() == 3 and  backbone_atoms[0].GetNeighbors()[0].GetAtomicNum() == 6 and backbone_atoms[0].GetNeighbors()[1].GetAtomicNum() == 8 and backbone_atoms[3].GetTotalDegree() == 3 and  backbone_atoms[3].GetNeighbors()[0].GetAtomicNum() == 6 and backbone_atoms[3].GetNeighbors()[1].GetAtomicNum() == 8):
          return False, "Molecule does not have CH2OH termini in linear chain"

    for i in range(1,len(backbone_atoms)-1):
         if not (backbone_atoms[i].GetTotalDegree() == 3 and  backbone_atoms[i].GetNeighbors()[0].GetAtomicNum() == 6 and backbone_atoms[i].GetNeighbors()[1].GetAtomicNum() == 6 and backbone_atoms[i].GetNeighbors()[2].GetAtomicNum() == 8 and backbone_atoms[i].GetNeighbors()[2].GetTotalDegree() == 1):
             return False, "Molecule does not have CHOH in linear chain"

    # Check that non-backbone carbon atoms have O neighbors
    for atom in atoms:
        if atom.GetAtomicNum() == 6 and atom not in backbone_atoms:
            has_oxygen_neighbor = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    has_oxygen_neighbor = True
                    break
            if not has_oxygen_neighbor:
                return False, "Non-backbone carbon atom does not have an oxygen neighbor"


    # Verify that there are only carbons, hydrogens, and oxygens
    for atom in atoms:
       if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Molecule contains atoms other than C, H, and O"

    return True, "Molecule is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH"