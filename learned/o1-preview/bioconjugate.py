"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize a set to store detected biological substructures
    bio_substructures = set()

    # Define SMARTS patterns for common biological substructures

    # Peptide bond (amide bond between amino acids)
    peptide_bond = Chem.MolFromSmarts("N[C;!$(C=O)]C(=O)")  # N-C-C(=O)
    if mol.HasSubstructMatch(peptide_bond):
        bio_substructures.add("peptide bond")

    # Sugar ring (five or six-membered ring with oxygen and multiple hydroxyl groups)
    sugar_ring = Chem.MolFromSmarts("C1OC([H])C([H])OC1")  # Simplified pattern for furanose ring
    if mol.HasSubstructMatch(sugar_ring):
        bio_substructures.add("sugar ring")

    # Nucleotide base (purine or pyrimidine ring systems)
    purine = Chem.MolFromSmarts("c1ncnc2ncnc12")  # Purine ring
    pyrimidine = Chem.MolFromSmarts("c1ncnc1")     # Pyrimidine ring
    if mol.HasSubstructMatch(purine) or mol.HasSubstructMatch(pyrimidine):
        bio_substructures.add("nucleotide base")

    # Phosphate group (common in nucleotides)
    phosphate_group = Chem.MolFromSmarts("P(=O)(O)O")  # Phosphate group
    if mol.HasSubstructMatch(phosphate_group):
        bio_substructures.add("phosphate group")

    # Long fatty acid chain (long aliphatic chain with terminal carboxylic acid)
    fatty_acid = Chem.MolFromSmarts("C(=O)OCCCCCCCCCC")  # 10+ carbon atoms
    if mol.HasSubstructMatch(fatty_acid):
        bio_substructures.add("fatty acid chain")

    # Glutathione substructure (tripeptide of glutamate, cysteine, and glycine)
    glutathione = Chem.MolFromSmarts("NCC(=O)N[C@@H](CS)C(=O)NCC(=O)O")  # Simplified glutathione pattern
    if mol.HasSubstructMatch(glutathione):
        bio_substructures.add("glutathione")

    # Coenzyme A substructure (contains ADP, pantetheine, and cysteamine)
    coa = Chem.MolFromSmarts("NC(=O)CCC(=O)NCCSC")  # Simplified CoA pattern
    if mol.HasSubstructMatch(coa):
        bio_substructures.add("coenzyme A")

    # Count the number of different biological substructures found
    num_bio_substructures = len(bio_substructures)

    # If at least two different biological substructures are found, classify as bioconjugate
    if num_bio_substructures >= 2:
        reason = f"Contains multiple biological substructures: {', '.join(bio_substructures)}"
        return True, reason
    else:
        return False, "Does not contain at least two different biological substructures"

# Example usage:
# smiles = "C(CCCC)[C@@H]([C@@H](\C=C\[C@H](C/C=C\C/C=C\CCCC(O)=O)O)SC[C@H](NC(CC[C@H](N)C(=O)O)=O)C(=O)NCC(=O)O)O"
# result = is_bioconjugate(smiles)
# print(result)