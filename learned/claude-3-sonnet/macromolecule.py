"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecules
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdqueries
from collections import Counter

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule has high molecular mass and consists of multiple repeating units.

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

    # Check molecular weight - typical macromolecules are >1000 Da
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a macromolecule"

    # Count atoms - macromolecules typically have many atoms
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 50:
        return False, f"Too few atoms ({num_atoms}) for a macromolecule"

    # Look for repeating structural patterns
    # Get all ring systems
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Count atom types and their frequencies
    atom_types = Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    
    # Count common functional groups
    patterns = {
        'amide': '[NX3][CX3](=[OX1])',
        'ester': '[OX2][CX3](=[OX1])',
        'ether': '[OX2][CX4]',
        'alcohol': '[OX2H]',
        'amine': '[NX3;H2,H1;!$(NC=O)]'
    }
    
    functional_groups = {}
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            matches = len(mol.GetSubstructMatches(pattern))
            functional_groups[name] = matches

    # Analyze repeating units based on functional groups
    total_functional_groups = sum(functional_groups.values())
    
    # Check for repeating units based on functional groups and rings
    has_repeating_units = (
        total_functional_groups >= 5 or  # Multiple functional groups
        num_rings >= 3 or  # Multiple ring systems
        (atom_types.get('C', 0) >= 20 and atom_types.get('O', 0) >= 10)  # Complex carbon skeleton
    )

    if not has_repeating_units:
        return False, "No significant repeating units found"

    # Calculate rotatable bonds - macromolecules typically have many
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, f"Too few rotatable bonds ({n_rotatable}) for a macromolecule"

    # Analyze complexity
    complexity_reasons = []
    if num_rings >= 3:
        complexity_reasons.append(f"{num_rings} ring systems")
    if total_functional_groups >= 5:
        complexity_reasons.append(f"{total_functional_groups} functional groups")
    for group, count in functional_groups.items():
        if count >= 3:
            complexity_reasons.append(f"{count} {group} groups")

    # Final classification
    if mol_wt > 1000 and has_repeating_units and len(complexity_reasons) >= 2:
        reason = f"Molecular weight: {mol_wt:.1f} Da, contains " + ", ".join(complexity_reasons)
        return True, reason
    
    return False, "Does not meet macromolecule criteria"