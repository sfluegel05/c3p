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

    # Determine which fragment is alcohol (RO-) and which is acid (RCOO-)
    def is_acid_fragment(frag):
        return any(atom.GetAtomicNum() == 8 and atom.GetDegree() == 2 for atom in frag.GetAtoms() if atom.GetIdx() != ester_matches[0][1])
    
    acid_frag, alcohol_frag = (fragments[0], fragments[1]) if is_acid_fragment(fragments[0]) else (fragments[1], fragments[0])

    # Check chain lengths using longest carbon chain
    def get_longest_chain(fragment):
        chains = []
        for atom in fragment.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1:  # Start from terminal carbon
                chain = rdMolDescriptors.CalcNumAtomStereoCenters(fragment, atom.GetIdx())
                chains.append(rdMolDescriptors.CalcLongestChain(fragment, atom.GetIdx()))
        return max(chains) if chains else 0

    alcohol_chain = get_longest_chain(alcohol_frag)
    acid_chain = get_longest_chain(acid_frag) - 1  # Subtract carboxyl carbon

    if alcohol_chain < 12 or acid_chain < 12:
        return False, f"Chains too short (alcohol: {alcohol_chain}, acid: {acid_chain})"

    # Check for other functional groups
    forbidden_patterns = [
        "[N,S,P]",  # No heteroatoms besides O in ester
        "[#6][F,Cl,Br,I]",  # No halogens
        "[*]=O",  # No other carbonyl groups
        "[*]=[*]=[*]",  # No allenes
        "[*]#[*]",  # No triple bonds
        "[*]O[C,S,N]",  # No ethers/thioesters/amides
    ]
    for pattern in forbidden_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains forbidden pattern: {pattern}"

    # Check molecular weight (typical wax esters >350 g/mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Contains single ester linking two long hydrocarbon chains"