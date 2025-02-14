"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def remove_inorganic_counterions(mol):
    """Removes common inorganic counterions from the molecule.
    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object
    Returns:
       rdkit.Chem.Mol: RDKit molecule object with counterions removed or None
    """
    if mol is None:
        return None
    
    # Create a copy of the molecule to avoid modifying the original
    new_mol = Chem.Mol(mol)
    
    # Identify and remove non-covalently bonded fragments
    for frag in Chem.GetMolFrags(new_mol, asMols=True):
        if all(atom.GetAtomicNum() not in [6,7,8,15,16] for atom in frag.GetAtoms()):
            new_mol = AllChem.DeleteSubstructs(new_mol, frag)

    return new_mol

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is likely a B vitamin based on its SMILES string.
    Checks for key structural features of B1, B2, B3, B5, B6, B7, B9 and B12 vitamins.
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if molecule is likely a B vitamin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Handle the special case of B12 directly from the name
    if 'B12' in smiles:
        return True, "Contains B12"

    # remove counterions and waters
    mol = remove_inorganic_counterions(mol)
    if mol is None:
        return False, "Could not remove counterions"

    # B1 (Thiamine) - Thiazole and pyrimidine rings
    thiamine_pattern = Chem.MolFromSmarts('c1[n,+]cscc1-CC-c1cnccnc1')
    if mol.HasSubstructMatch(thiamine_pattern):
       return True, "Contains thiamine (B1) substructure."


    # B2 (Riboflavin) - Isoalloxazine ring system
    riboflavin_pattern = Chem.MolFromSmarts('C1=CC=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1')
    if mol.HasSubstructMatch(riboflavin_pattern):
        return True, "Contains riboflavin (B2) substructure."

    # B3 (Niacin/Nicotinamide) - Pyridine with carboxyl/amide
    nicotinamide_pattern = Chem.MolFromSmarts('c1ccncc1C(=O)[N]')
    nicotinic_acid_pattern = Chem.MolFromSmarts('c1ccncc1C(=O)O')
    if mol.HasSubstructMatch(nicotinamide_pattern) or mol.HasSubstructMatch(nicotinic_acid_pattern):
          return True, "Contains niacin/nicotinamide (B3) substructure."

    # B5 (Pantothenic acid) - Pantoic acid and beta-alanine
    pantothenic_pattern = Chem.MolFromSmarts('CC(C)(CO)[C@H](O)C(=O)NCCC(O)=O')
    if mol.HasSubstructMatch(pantothenic_pattern):
       return True, "Contains pantothenic acid (B5) substructure."
    
    # B6 (Pyridoxine, Pyridoxal, Pyridoxamine) - Pyridine with OH, CH2OH, CHO, and/or CH2NH2 groups
    pyridoxine_core_pattern = Chem.MolFromSmarts('c1c(O)cnc(C)c1[C]')
    if mol.HasSubstructMatch(pyridoxine_core_pattern):
       return True, "Contains pyridoxine/pyridoxal/pyridoxamine (B6) substructure."

    # B7 (Biotin) - Tetrahydrothiophene and imidazole
    biotin_pattern = Chem.MolFromSmarts('C1S[CH]2[CH]([CH2][C](N2)=O)NC(=O)N[C@H]3CCCC[C@H]3')
    if mol.HasSubstructMatch(biotin_pattern):
        return True, "Contains biotin (B7) substructure."

    # B9 (Folic acid/folate) - Pteridine with PABA and glutamic acid
    folic_acid_core_pattern = Chem.MolFromSmarts('c1nc2nc(N)[nH]c(=O)c2n1-CC-c1ccc(N)cc1-C(=O)NC(CCC(=O)O)C(=O)O')
    pteridine_pattern = Chem.MolFromSmarts('c1nc2nc(N)[nH]c(=O)c2n1')
    if mol.HasSubstructMatch(folic_acid_core_pattern) or mol.HasSubstructMatch(pteridine_pattern):
        return True, "Contains folic acid (B9) substructure."


    return False, "Not recognized as a B vitamin based on substructure search."