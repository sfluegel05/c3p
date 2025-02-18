"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:100157 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is an ester resulting from a fatty acid and a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one ester group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Split molecule into two parts at the ester bond
    ester_bond = mol.GetBondBetweenAtoms(ester_matches[0][0], ester_matches[0][1])
    if not ester_bond:
        return False, "Invalid ester bond structure"
    
    # Break the ester bond to separate alcohol and acid parts
    fragmented_mols = Chem.FragmentOnBonds(mol, [ester_bond.GetIdx()], addDummies=False)
    fragments = Chem.GetMolFrags(fragmented_mols, asMols=True)
    if len(fragments) != 2:
        return False, "Ester does not split molecule into two parts"

    # Determine which fragment is acid (RCOO-) and which is alcohol (RO-)
    def is_acid_fragment(frag):
        # Acid fragment has the carbonyl oxygen
        return any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0 for atom in frag.GetAtoms())
    
    acid_frag, alcohol_frag = (fragments[0], fragments[1]) if is_acid_fragment(fragments[0]) else (fragments[1], fragments[0])

    # Calculate longest carbon chain for each fragment
    acid_chain = rdMolDescriptors.CalcLongestChain(acid_frag)
    alcohol_chain = rdMolDescriptors.CalcLongestChain(alcohol_frag)

    # Adjust acid chain length (subtract 1 for the carbonyl carbon in the ester)
    acid_chain -= 1
    if acid_chain < 8:
        return False, f"Acid chain too short ({acid_chain} carbons)"
    if alcohol_chain < 8:
        return False, f"Alcohol chain too short ({alcohol_chain} carbons)"

    # Check for exactly two oxygen atoms (ester group only)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, f"Expected exactly 2 oxygen atoms, found {o_count}"

    # Check for absence of other heteroatoms (N, S, P, halogens)
    heteroatoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 8, 1]]
    if heteroatoms:
        return False, f"Contains forbidden heteroatoms: {heteroatoms}"

    # Check molecular weight (typical wax esters are >350 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Single ester of long-chain fatty acid and fatty alcohol"