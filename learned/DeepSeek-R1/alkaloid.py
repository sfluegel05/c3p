"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:22333 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, MolFromSmarts
from rdkit.Chem.rdMolDescriptors import CalcNumAmideBonds

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are nitrogen-containing compounds with heterocyclic structures,
    excluding amino acids, peptides, proteins, nucleotides, and exocyclic amines.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Basic nitrogen check
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        return False, "No nitrogen atoms"

    # Exclude amino acids (alpha-amine + alpha-carboxylic acid)
    amino_acid_pattern = MolFromSmarts("[NH2]-[CX4]-C(=O)[OH]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Amino acid structure detected"

    # Exclude peptides (amide bonds)
    if CalcNumAmideBonds(mol) > 0:
        return False, "Peptide bonds present"

    # Exclude nucleotides (phosphate + ribose-like structure)
    phosphate = MolFromSmarts("[P](=O)(O)(O)O")
    ribose = MolFromSmarts("C1OCC(O)C1")
    if mol.HasSubstructMatch(phosphate) and mol.HasSubstructMatch(ribose):
        return False, "Nucleotide detected"

    # Check for heterocyclic nitrogen (any nitrogen in a ring)
    ring_info = mol.GetRingInfo()
    heterocyclic_n = any(ring_info.NumAtomRings(n.GetIdx()) > 0 for n in nitrogens)
    
    if heterocyclic_n:
        # Check for excluded heterocycles (pyrimidine/purine in nucleotides)
        pyrimidine = MolFromSmarts("n1cncc1")  # Simplified pyrimidine pattern
        purine = MolFromSmarts("n1c2ncnc2cn1")
        if not (mol.HasSubstructMatch(pyrimidine) or mol.HasSubstructMatch(purine)):
            return True, "Heterocyclic nitrogen present"
    
    # Exclude simple exocyclic amines with low molecular weight
    amine_pattern = MolFromSmarts("[NH2,NH1]")
    if mol.HasSubstructMatch(amine_pattern):
        mol_wt = Descriptors.ExactMolWt(mol)
        if mol_wt < 250:  # Adjusted threshold
            return False, "Low molecular weight amine"

    # Check for complex ring systems with nitrogen (even if not in ring)
    if ring_info.NumRings() > 1 and len(nitrogens) >= 1:
        return True, "Complex multi-ring system with nitrogen"

    return False, "Does not meet alkaloid criteria"