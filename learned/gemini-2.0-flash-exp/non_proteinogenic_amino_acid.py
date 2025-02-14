"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is an amino acid not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Basic amino acid structure check: N-C-C(=O)-O and a sidechain
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4]([H])[CX3](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Missing basic amino acid structure"

    # 2. Exclude standard proteinogenic amino acids (side chains)

    # Common side chains (simple alkyl, aryl, and hydroxylated)
    simple_side_chain_patterns = [
        Chem.MolFromSmarts("[CX4H3]"), # Methyl
        Chem.MolFromSmarts("[CX4H2][CX4H3]"), # Ethyl
        Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H3]"), # Propyl
        Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H3]"), # Butyl
        Chem.MolFromSmarts("[CX4H2]1[CX4H2][CX4H2][CX4H2]1"), # Cyclobutyl
        Chem.MolFromSmarts("[CX4H2]1[CX4H2][CX4H2]1"), # Cyclopropyl
        Chem.MolFromSmarts("c1ccccc1"), # Phenyl
        Chem.MolFromSmarts("c1ccncc1"), # Pyridyl
        Chem.MolFromSmarts("[CX4H2][OX2H1]"), # Hydroxymethyl
        Chem.MolFromSmarts("[CX4H2][CX4H2][OX2H1]"), # Hydroxyethyl
        Chem.MolFromSmarts("CC(C)"), # Isopropyl
        Chem.MolFromSmarts("CC(C)C"), # Isobutyl
        Chem.MolFromSmarts("CC(C)(C)"), # tert-butyl
        Chem.MolFromSmarts("C[C@@H](C)"), # chiral version of isopropyl
        Chem.MolFromSmarts("c1cc[nH]c1"), # Imidazole
        Chem.MolFromSmarts("C[S]"), # methylthio
        Chem.MolFromSmarts("SC"), # methylthio, alternative notation
        Chem.MolFromSmarts("CC(=O)N"), # Acetamide, present in asparagine and glutamine

    ]

    # Check if any simple side chain is directly attached to the alpha carbon of the amino acid.
    # If one of them is attached, then assume it is not a non-canonical amino acid
    alpha_carbon_with_sidechain_pattern = Chem.MolFromSmarts("[NX3][CX4]([H])([!H])[CX3](=[OX1])[OX2]")
    alpha_carbon_match = mol.GetSubstructMatch(alpha_carbon_with_sidechain_pattern)
    if alpha_carbon_match:
      alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon_match[1]) # Get alpha carbon
      for pattern in simple_side_chain_patterns:
        for atom in mol.GetAtoms():
          if atom.GetIdx() != alpha_carbon_atom.GetIdx() and atom.GetIdx() in [neighbor.GetIdx() for neighbor in alpha_carbon_atom.GetNeighbors()]:
            submol = Chem.PathToSubmol(mol, [atom.GetIdx()] + [neighbor.GetIdx() for neighbor in atom.GetNeighbors()])
            if submol.HasSubstructMatch(pattern):
              return False, "Has a common proteinogenic amino acid side chain"


    # 3. Check for Modifications
    modification_patterns = [
        Chem.MolFromSmarts("[Cl]"),  # Chlorine
        Chem.MolFromSmarts("[Br]"),  # Bromine
        Chem.MolFromSmarts("[F]"),  # Fluorine
        Chem.MolFromSmarts("[IX1]"), # Iodine
        Chem.MolFromSmarts("[OX2H]"), # Extra hydroxyl group
        Chem.MolFromSmarts("O[P](=O)(O)O"), # Phosphate group
        Chem.MolFromSmarts("[CX3](=[OX1])[NX2]"), # Amide group
        Chem.MolFromSmarts("[NX2]O"), # Hydroxylamine
        Chem.MolFromSmarts("[B](O)O"), # Boronic acid
        Chem.MolFromSmarts("[CX3]#[NX1]")  # Cyano group
    ]

    has_modification = False
    for pattern in modification_patterns:
      if mol.HasSubstructMatch(pattern):
        has_modification = True
        break
    
    # 4. Combine Checks
    if has_modification:
      return True, "Contains modifications typical of non-proteinogenic amino acids"
    else:
        return True, "Does not have common side chain; but no specific modifications identified"