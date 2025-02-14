"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: monocarboxylic acid anion
"""
from rdkit import Chem

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove counter ions and separate molecules
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    largest_mol = max(fragments, default=mol, key=lambda m: m.GetNumAtoms())

    # Use the largest fragment for classification
    mol = largest_mol

    # Define SMARTS patterns
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[O-,OH]")  # Protonated or deprotonated carboxyl groups
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")   # Deprotonated carboxylate group
    charged_groups_pattern = Chem.MolFromSmarts("[+,-]")        # Any atom with a formal charge
    metal_pattern = Chem.MolFromSmarts("[Na,K,Ca,Mg,Zn,Cu,Fe,Co,Ni,Al,Li,Ag,Hg,Pb,Tl]")  # Metals

    # Check for metals (exclude metal salts)
    if mol.HasSubstructMatch(metal_pattern):
        return False, "Molecule contains metal ions"

    # Find all carboxyl groups
    carboxyl_groups = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl_groups = len(carboxyl_groups)
    
    if num_carboxyl_groups != 1:
        return False, f"Found {num_carboxyl_groups} carboxyl groups, expected exactly 1"

    # Find deprotonated carboxylate groups
    carboxylate_groups = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate_groups = len(carboxylate_groups)
    
    if num_carboxylate_groups != 1:
        return False, "Carboxyl group is not deprotonated"

    # Check for other charged groups
    charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]
    # Exclude charges on carboxylate oxygen
    carboxylate_oxygens = [mol.GetAtomWithIdx(match[1]) for match in carboxylate_groups]
    other_charged_atoms = [atom for atom in charged_atoms if atom not in carboxylate_oxygens]
    
    if len(other_charged_atoms) > 0:
        return False, "Molecule contains other charged groups"

    # Calculate net charge
    net_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if net_charge != -1:
        return False, f"Net charge is {net_charge}, expected -1"

    # If all conditions are met, classify as monocarboxylic acid anion
    return True, "Molecule is a monocarboxylic acid anion with exactly one deprotonated carboxyl group"