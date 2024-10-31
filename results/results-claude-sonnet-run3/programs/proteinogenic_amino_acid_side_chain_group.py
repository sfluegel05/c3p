from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid_side_chain_group(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid side chain group.
    
    Args:
        smiles (str): SMILES string of the molecule, with * marking attachment point
        
    Returns:
        bool: True if molecule is a proteinogenic amino acid side chain group
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
            
        # Check if molecule contains exactly one attachment point (*)
        attachment_points = len([a for a in mol.GetAtoms() if a.GetSymbol() == '*'])
        if attachment_points != 1:
            return False, f"Molecule must contain exactly one attachment point (*), found {attachment_points}"

        # Get attachment point atom
        attachment_atom = next(a for a in mol.GetAtoms() if a.GetSymbol() == '*')
        
        # Check that attachment point has exactly one neighbor
        if len(attachment_atom.GetNeighbors()) != 1:
            return False, "Attachment point must have exactly one neighbor"
            
        # Check that neighbor of attachment point is carbon
        neighbor = attachment_atom.GetNeighbors()[0]
        if neighbor.GetSymbol() != 'C':
            return False, "Attachment point must be connected to carbon"

        # Check for disallowed features that would indicate non-proteinogenic side chains
        # No triple bonds allowed
        if any(bond.GetBondType() == Chem.BondType.TRIPLE for bond in mol.GetBonds()):
            return False, "Triple bonds not allowed in proteinogenic side chains"
            
        # No halogen atoms allowed
        if any(atom.GetSymbol() in ['F', 'Cl', 'Br', 'I'] for atom in mol.GetAtoms()):
            return False, "Halogens not allowed in proteinogenic side chains"
            
        # No double bonds to attachment point
        for bond in attachment_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                return False, "Double bonds to attachment point not allowed"

        # List of allowed atoms in proteinogenic side chains
        allowed_atoms = {'C', 'H', 'N', 'O', 'S', 'Se', '*'}
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in allowed_atoms:
                return False, f"Atom {atom.GetSymbol()} not allowed in proteinogenic side chains"

        # If we get here, molecule has passed all checks
        return True, "Valid proteinogenic amino acid side chain group"
        
    except Exception as e:
        return None, f"Error processing molecule: {str(e)}"
# Pr=0.8260869565217391
# Recall=1.0