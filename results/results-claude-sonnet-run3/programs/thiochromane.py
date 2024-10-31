from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_thiochromane(smiles: str):
    """
    Determines if a molecule is a thiochromane (benzene ring fused to a thiopyran).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiochromane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for presence of sulfur
    if not any(atom.GetSymbol() == 'S' for atom in mol.GetAtoms()):
        return False, "No sulfur atom found"
        
    # Get ring info
    rings = mol.GetRingInfo()
    if not rings.NumRings():
        return False, "No rings found"

    # Find 6-membered aromatic rings (benzene)
    benzene_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() and atom.GetSymbol() == 'C' for atom in atoms):
                benzene_rings.append(ring)

    if not benzene_rings:
        return False, "No benzene ring found"

    # Find rings containing sulfur
    sulfur_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetSymbol() == 'S' for atom in atoms):
            sulfur_rings.append(ring)

    if not sulfur_rings:
        return False, "No ring containing sulfur found"

    # Check if benzene and sulfur rings share atoms (are fused)
    for benzene_ring in benzene_rings:
        for sulfur_ring in sulfur_rings:
            shared_atoms = set(benzene_ring).intersection(set(sulfur_ring))
            if len(shared_atoms) == 2:  # Two shared atoms indicate fusion
                # Verify sulfur ring is 6-membered (thiopyran) or 5-membered (thiophene)
                if len(sulfur_ring) in [5,6]:
                    sulfur_atom = next(mol.GetAtomWithIdx(i) for i in sulfur_ring 
                                    if mol.GetAtomWithIdx(i).GetSymbol() == 'S')
                    
                    if sulfur_atom.GetDegree() == 2:  # Basic thiopyran/thiophene
                        return True, "Contains benzene fused to thiopyran/thiophene ring"
                    elif sulfur_atom.GetDegree() > 2:  # Oxidized forms also count
                        return True, "Contains benzene fused to oxidized thiopyran/thiophene ring"

    return False, "No fused benzene-thiopyran system found"
# Pr=1.0
# Recall=1.0