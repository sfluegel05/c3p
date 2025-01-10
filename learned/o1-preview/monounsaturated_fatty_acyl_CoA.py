"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: monounsaturated fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    '''
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is a fatty acyl-CoA in which the fatty acyl chain contains exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    '''

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'Invalid SMILES string'

    # Check for thioester linkage [C(=O)S]
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, 'No thioester linkage found'

    # Check for adenine ring (as part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, 'No adenine moiety found (not a CoA derivative)'

    # Assume the first thioester is the fatty acyl linkage
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]  # Carbonyl carbon index
    sulfur_idx = thioester_match[2]      # Sulfur atom index

    # Break the bond between carbonyl carbon and sulfur to isolate the acyl chain
    emol = Chem.EditableMol(mol)
    emol.RemoveBond(carbonyl_c_idx, sulfur_idx)
    fragmented_mol = emol.GetMol()

    # Get the fragments - acyl chain and CoA moiety
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)

    # Identify the acyl chain fragment (should contain the carbonyl carbon)
    acyl_chain = None
    coa_moiety = None
    for frag in fragments:
        atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
        if carbonyl_c_idx in atom_indices:
            acyl_chain = frag
        elif sulfur_idx in atom_indices:
            coa_moiety = frag

    if acyl_chain is None:
        return False, 'Could not isolate acyl chain'

    # Check that the acyl chain is linear and aliphatic
    if acyl_chain.GetRingInfo().NumRings() > 0:
        return False, 'Acyl chain contains rings (not linear)'

    # Count number of carbon atoms in the acyl chain (excluding carbonyl carbon)
    carbon_atoms = [atom for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 4:
        return False, 'Acyl chain too short to be a fatty acid'

    # Count carbon-carbon double bonds in the acyl chain
    num_double_bonds = 0
    for bond in acyl_chain.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                num_double_bonds += 1

    if num_double_bonds != 1:
        return False, f'Acyl chain contains {num_double_bonds} carbon-carbon double bonds (expected exactly 1)'

    return True, 'Molecule is a monounsaturated fatty acyl-CoA (contains exactly one carbon-carbon double bond in the acyl chain)'

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.'
    }
}