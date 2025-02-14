"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:27270 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic amino compound with an amino group and a hydroxy group
    attached to the same carbon atom, which is not a hydrogen atom. It is an intermediate
    in the formation of imines from aldehydes or ketones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for atoms with both amino (-NH2/-NH-) and hydroxy (-OH) groups attached to the same carbon
    hemiaminal_pattern = Chem.MolFromSmarts("[NH2,NH1][C;H1]([OH1])[!#1]")
    hemiaminal_atoms = mol.GetSubstructMatches(hemiaminal_pattern)
    
    if not hemiaminal_atoms:
        return False, "No atoms with amino and hydroxy groups attached to the same carbon atom"
    
    # Check if the molecule contains a carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)[!#1]")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    # Assign stereochemistry and identify unique molecules
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    inchi_key = AllChem.AssignBondEndingUniqueCopyOfInchiKey(mol)
    
    if has_carbonyl:
        return True, f"Contains amino and hydroxy groups attached to the same carbon atom, and a carbonyl group (InChI Key: {inchi_key})"
    else:
        return False, f"Does not contain a carbonyl group (InChI Key: {inchi_key})"