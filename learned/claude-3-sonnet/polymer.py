"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: Polymers (CHEBI:24702)

A polymer is a mixture, which is composed of macromolecules of different kinds
and which may be differentiated by composition, length, degree of branching, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - polymers typically have high MW
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a polymer"

    # Look for repeating units/patterns in the SMILES
    smiles_parts = smiles.split('.')
    if len(smiles_parts) > 1:
        # Multiple components, could be a mixture of macromolecules
        repeating_parts = [part for part in smiles_parts if smiles.count(part) > 1]
        if repeating_parts:
            return True, "Contains repeating structural units, could be a polymer mixture"

    # Check for long carbon chains or rings (>6 atoms)
    rings = mol.GetRingInfo().AtomRings()
    large_rings = [len(ring) for ring in rings if len(ring) > 6]
    chains = [len(chain) for chain in Chem.FindAllPathsOfLengthN(mol, 6, useBonds=False)]
    if large_rings or chains:
        return True, "Contains large rings or long carbon chains, potentially a polymer"

    # No definitive evidence of polymeric structure
    return False, "No clear indications of a polymeric structure"