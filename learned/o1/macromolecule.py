"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS

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

    # Generate a list of fragments (attempting to find repeating units)
    # Here we use the RDKit function to enumerate possible ring systems or substructures
    # This is a simplistic approach and may not capture all repeating units
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) <= 1:
        # If the molecule is a single fragment, we attempt to decompose it
        # into simpler units using BRICS decomposition
        brics_fragments = Chem.BRICS.BRICSDecompose(mol)
        if len(brics_fragments) <= 1:
            return False, "No repeating units detected in the molecule"

        # Check for repeating fragments
        frag_counts = {}
        for frag_smiles in brics_fragments:
            frag_counts[frag_smiles] = frag_counts.get(frag_smiles, 0) + 1

        repeating_units = {smiles: count for smiles, count in frag_counts.items() if count > 1}
        if not repeating_units:
            return False, "No repeating units detected in the molecule"
    else:
        # Multiple disconnected fragments detected
        # Check for repeating fragments
        frag_smiles_list = [Chem.MolToSmiles(frag) for frag in frags]
        frag_counts = {}
        for frag_smiles in frag_smiles_list:
            frag_counts[frag_smiles] = frag_counts.get(frag_smiles, 0) + 1

        repeating_units = {smiles: count for smiles, count in frag_counts.items() if count > 1}
        if not repeating_units:
            return False, "No repeating units detected in the molecule"

    # If we reach here, the molecule is large and has repeating units
    return True, f"Molecule is a macromolecule with molecular weight {mol_wt:.2f} Da and repeating units detected"