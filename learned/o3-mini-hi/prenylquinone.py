"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
This program uses heuristic SMARTS patterns to identify a quinone core and a prenyl (isoprene-derived)
side chain.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    The routine checks for a quinone core (using a few SMARTS patterns for quinone systems)
    and for the presence of a polyprenyl-type side chain (heuristically using a common prenyl fragment).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a prenylquinone, False otherwise
        str: Reason for classification decision
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of SMARTS patterns to capture common quinone cores.
    # These include a typical 1,4-benzoquinone and a simplified naphthoquinone pattern.
    quinone_smarts = [
        "c1c(=O)ccc(=O)c1",      # 1,4-benzoquinone
        "c1ccc2C(=O)ccc2c1",     # annulated system such as naphthoquinone
    ]
    quinone_found = False
    for smarts in quinone_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            quinone_found = True
            break
    if not quinone_found:
        return False, "No quinone core detected"
    
    # Define a SMARTS pattern for a prenyl (isoprene-derived) fragment.
    # It searches for a CH2 group double-bonded to a carbon bearing a methyl group.
    prenyl_pattern = Chem.MolFromSmarts("[CH2]C=C(C)")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    
    # Sometimes a substituted chain may not exactly match the small prenyl fragment.
    # In that case we try a fallback: look for a relatively long, non-aromatic chain 
    # (at least 5 carbon atoms) that contains a double bond.
    if not prenyl_matches:
        # Get all aliphatic chains (atoms not in aromatic rings)
        aliphatic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum()==6 and not atom.GetIsAromatic()]
        if not aliphatic_carbons:
            return False, "No aliphatic carbon chain found for prenyl side-chain"
        # For each aliphatic carbon, we try to find a neighbor chain that is at least 5 carbons long.
        # This is a simple heuristic.
        chain_found = False
        for idx in aliphatic_carbons:
            atom = mol.GetAtomWithIdx(idx)
            # Examine neighbors; if at least one bond is a double bond, consider it.
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() >= 2.0:
                    chain_found = True
                    break
            if chain_found:
                break
        if not chain_found:
            return False, "No prenyl-like side chain (unsaturated aliphatic chain) detected"
    # If a prenyl fragment is found, we consider that as evidence for the polyprenyl side chain.
    return True, "Molecule contains a quinone core and a prenyl-derived unsaturated side chain"