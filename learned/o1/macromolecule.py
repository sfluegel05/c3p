"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import BRICS
from rdkit.Chem import rdFingerprintGenerator

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule of high relative molecular mass,
    the structure of which essentially comprises the multiple repetition of units
    derived from molecules of low relative molecular mass.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight is {mol_wt:.2f} Da, which is below the macromolecule threshold"

    # Generate molecular fingerprint
    fp_generator = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fp = fp_generator.GetFingerprint(mol)

    # Count subgraph frequencies
    substructures = {}
    for atom in mol.GetAtoms():
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=1, atomIdx=atom.GetIdx())
        submol = Chem.PathToSubmol(mol, env)
        smi = Chem.MolToSmiles(submol, canonical=True)
        if smi:
            substructures[smi] = substructures.get(smi, 0) + 1

    # Identify repeating units
    repeating_units = {smi: count for smi, count in substructures.items() if count > 1}

    if not repeating_units:
        return False, "No repeating units detected in the molecule"

    # If we reach here, the molecule is large and has repeating units
    return True, f"Molecule is a macromolecule with molecular weight {mol_wt:.2f} Da and repeating units detected"