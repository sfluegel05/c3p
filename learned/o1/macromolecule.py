"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from collections import Counter

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
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Calculate molecular weight
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 500:
            return False, f"Molecular weight is {mol_wt:.2f} Da, which is below typical macromolecule mass"

        # Improve bond patterns
        peptide_bond = Chem.MolFromSmarts("N[C;!$(C=O)][C;!$(C=O)](=O)")
        glycosidic_bond = Chem.MolFromSmarts("[C;!H0;!$(C=O)][O][C;!H0;!$(C=O)]")  # Exclude carbonyl carbons
        ester_bond = Chem.MolFromSmarts("C(=O)O[C;!$(C=O)]")
        amide_bond = Chem.MolFromSmarts("C(=O)N[C;!$(C=O)]")

        bond_patterns = [peptide_bond, glycosidic_bond, ester_bond, amide_bond]
        pattern_names = ['peptide bonds', 'glycosidic bonds', 'ester bonds', 'amide bonds']
        repeating_units_detected = []

        # Check for repeating substructures
        frag_mols = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        frag_smiles = [Chem.MolToSmiles(frag, isomericSmiles=True) for frag in frag_mols]
        frag_counts = Counter(frag_smiles)

        # Find fragments that occur multiple times
        repeating_frags = [frag for frag, count in frag_counts.items() if count > 1]

        # Look for repeating bond patterns
        total_repeats = 0
        for pattern, name in zip(bond_patterns, pattern_names):
            matches = mol.GetSubstructMatches(pattern)
            if len(matches) >= 2:  # Threshold for considering repeating units
                repeating_units_detected.append(name)
                total_repeats += len(matches)

        if not repeating_units_detected and not repeating_frags:
            return False, "No repeating units detected in the molecule"

        # Combine reasons
        reasons = f"Molecule is a macromolecule with molecular weight {mol_wt:.2f} Da"
        if repeating_units_detected:
            reasons += f" and contains repeating units: {', '.join(repeating_units_detected)}"
        if repeating_frags:
            reasons += f" and has repeating substructures: {', '.join(set(repeating_frags))}"

        return True, reasons

    except Exception as e:
        return False, f"An error occurred during processing: {str(e)}"