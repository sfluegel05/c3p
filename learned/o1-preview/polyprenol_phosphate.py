"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation 
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phosphate group pattern (including mono-, di-, tri-phosphates)
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(O)(O)=O"),            # Monophosphate
        Chem.MolFromSmarts("OP(O)(=O)OP(O)(O)=O"),   # Diphosphate
        Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)OP(O)(O)=O")  # Triphosphate
    ]

    # Search for phosphate groups
    phosphate_matches = []
    for pattern in phosphate_patterns:
        matches = mol.GetSubstructMatches(pattern)
        phosphate_matches.extend(matches)

    if not phosphate_matches:
        return False, "No phosphate group found"

    # Define a general isoprene unit pattern: C=C-C(-C)
    isoprene_smarts = "[CH2]=[CH][CH2][CH2]"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)

    # Function to recursively traverse the molecule from a starting atom
    def traverse(atom, visited, chain_atoms):
        visited.add(atom.GetIdx())
        chain_atoms.append(atom)
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetIdx() not in visited:
                # Continue traversal through carbon, oxygen, or phosphorus atoms
                if neighbor.GetAtomicNum() in [6, 8, 15]:
                    traverse(neighbor, visited, chain_atoms)

    # Trace from each phosphate group
    for phosphate_match in phosphate_matches:
        phosphate_atoms = [mol.GetAtomWithIdx(idx) for idx in phosphate_match]
        # Start traversal from the oxygen atom connected to the phosphate group
        for atom in phosphate_atoms:
            if atom.GetAtomicNum() == 8:
                visited = set()
                chain_atoms = []
                traverse(atom, visited, chain_atoms)
                # Extract carbon chain atoms
                carbon_chain_atoms = [a for a in chain_atoms if a.GetAtomicNum() == 6]
                # Check if there is a long carbon chain
                if len(carbon_chain_atoms) >= 20:
                    # Check for isoprene units within the chain
                    submol = Chem.PathToSubmol(mol, [a.GetIdx() for a in carbon_chain_atoms])
                    isoprene_matches = submol.GetSubstructMatches(isoprene_pattern)
                    num_isoprene_units = len(isoprene_matches)
                    if num_isoprene_units >= 4:
                        return True, f"Molecule is a polyprenol phosphate with {num_isoprene_units} isoprene units"
                    else:
                        # Alternative: Estimate isoprene units by counting double bonds
                        num_double_bonds = sum(1 for bond in submol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
                        if num_double_bonds >= 4:
                            return True, f"Molecule is a polyprenol phosphate with {num_double_bonds} double bonds (estimated isoprene units)"
                        else:
                            continue
    return False, "No polyprenol chain connected to phosphate group with sufficient isoprene units"