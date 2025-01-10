"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a large macrocyclic lactone ring.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify if the molecule has a macrolide ring.
    ring_info = mol.GetRingInfo()
    ring_bonds = ring_info.BondRings()
    solvent_heavy_fragments = list(Chem.GetMolFrags(mol, sanitizeFrags=False, asMols=True))

    # Identifying candidate macrolide rings
    for ring in ring_bonds:
        if len(ring) >= 12 and len(ring) <= 16:
            ring_atoms = [mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx() for bond_idx in ring] + \
                         [mol.GetBondWithIdx(ring[-1]).GetEndAtomIdx()]
            submol = rdmolops.FragmentOnBonds(mol, bondsToBreak=[], atomsToKeep=ring_atoms)

            # Look for ester linkage within the ring
            submol_smarts = Chem.MolToSmiles(submol)
            if "C(=O)O" in submol_smarts:
                return True, "Contains macrocyclic lactone ring with ester linkage"

    return False, "No characteristic macrolide macrocyclic lactone ring found"