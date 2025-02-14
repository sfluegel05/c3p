"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea is a urea group where one nitrogen is bonded to a sulfonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for urea
    urea_pattern = Chem.MolFromSmarts("NC(=O)N")
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "Does not contain a urea group"

    # Define the SMARTS pattern for sulfonyl attached to nitrogen
    sulfonyl_pattern = Chem.MolFromSmarts("NS(=O)(=O)")
    #Get all nitrogens connected to a nitrogen in the urea group
    urea_nitrogen_matches = mol.GetSubstructMatches(urea_pattern)
    
    for match in urea_nitrogen_matches:
        urea_nitrogens = [match[0], match[2]]
        for urea_nitrogen in urea_nitrogens:
            for atom in mol.GetAtomWithIdx(urea_nitrogen).GetNeighbors():
                if atom.GetAtomicNum() == 7:
                    temp_mol = Chem.RWMol(mol)
                    bond = temp_mol.GetBondBetweenAtoms(urea_nitrogen, atom.GetIdx())
                    bond_order = bond.GetBondType()
                    temp_mol.RemoveBond(urea_nitrogen, atom.GetIdx())
                    
                    temp_smiles = Chem.MolToSmiles(temp_mol)
                    
                    temp_mol2 = Chem.MolFromSmiles(temp_smiles)
                    if temp_mol2:
                        sub_mol = Chem.MolFromSmarts(f"[{atom.GetIdx()}]NS(=O)(=O)")
                        if temp_mol2.HasSubstructMatch(sub_mol):
                            return True, "Contains N-sulfonylurea group"
                    
                    bond = temp_mol.AddBond(urea_nitrogen, atom.GetIdx(), bond_order)

    return False, "Does not contain N-sulfonylurea group"