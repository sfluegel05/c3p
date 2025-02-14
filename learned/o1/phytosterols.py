"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:18374 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    A phytosterol is a plant sterol similar to cholesterol, varying only in carbon side chains
    and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general sterol core SMARTS pattern (steroid nucleus)
    # Four fused rings: three six-membered rings and one five-membered ring
    sterol_core_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -[#6]2-[#6]=[#6]-[#6]-[#6]-[#6]-2
    -[#6]3-[#6]-[#6]-[#6]-[#6]-3
    -[#6]4-[#6]-[#6]-[#6]-[#6]-4
    """

    sterol_core_smarts = sterol_core_smarts.replace("\n", "").replace(" ", "")
    sterol_core = Chem.MolFromSmarts(sterol_core_smarts)
    if sterol_core is None:
        return False, "Error in sterol core SMARTS pattern"

    # Check if molecule contains sterol core
    if not mol.HasSubstructMatch(sterol_core):
        return False, "Steroid nucleus not found"

    # Check for hydroxyl group at position 3
    # Position 3 is adjacent to ring A of the steroid nucleus
    hydroxyl_smarts = "[#6R1]-[#6R1]-[#6R1]([#8H])"  # Simplified pattern for C-C-C(OH)
    hydroxyl_group = Chem.MolFromSmarts(hydroxyl_smarts)
    if not mol.HasSubstructMatch(hydroxyl_group):
        return False, "No hydroxyl group at position 3"

    # Optional: Check that variations are only in side chains and double bonds
    # Generate the steroid core of the molecule
    mol_core = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[!R]"))
    mol_core = Chem.RemoveHs(mol_core)

    # Generate the steroid core of cholesterol for comparison
    cholesterol_smiles = 'C[C@H](CCCC(C)C)C1CCC2C1(C)CC=C3C2CCC4[C@@H](O)CC[C@]34C'
    cholesterol_mol = Chem.MolFromSmiles(cholesterol_smiles)
    cholesterol_core = Chem.DeleteSubstructs(cholesterol_mol, Chem.MolFromSmarts("[!R]"))
    cholesterol_core = Chem.RemoveHs(cholesterol_core)

    # Compare cores using RDKit's Isomorphism
    from rdkit.Chem import rdFMCS
    res = rdFMCS.FindMCS([mol_core, cholesterol_core], ringMatchesRingOnly=True,
                         completeRingsOnly=True, matchChiralTag=False)

    if res.numAtoms < 17:  # Steroid nucleus has 17 carbon atoms
        return False, "Core structure differs significantly from cholesterol"

    # If all checks pass, classify as phytosterol
    return True, "Molecule contains steroid nucleus and matches phytosterol criteria"