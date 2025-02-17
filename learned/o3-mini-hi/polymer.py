"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is a mixture composed of macromolecules of different kinds,
which may be differentiated by composition, length, degree of branching etc..
This program implements a simple heuristic for classifying a SMILES string as polymer.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer (i.e. a mixture
    of macromolecules of different kinds). This heuristic splits the molecule into
    disconnected fragments, ignores very small fragments (such as ions) and then
    checks if there are at least two different large fragments.
    
    Args:
        smiles (str): SMILES string of the substance
        
    Returns:
        bool: True if classified as a polymer, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split the molecule into fragments (each fragment is a disconnected component)
    frags = Chem.GetMolFrags(mol, asMols=True)
    
    # We ignore very small fragments (e.g. common counter ions) by filtering on molecular weight.
    # Here, a threshold of 100 Da is used (this can be tuned).
    large_frags = [frag for frag in frags if rdMolDescriptors.CalcExactMolWt(frag) >= 100]
    
    # If we have two or more large fragments, then check if they are chemically different.
    # A polymer (by our definition) is a mixture of macromolecules of different kinds.
    if len(large_frags) >= 2:
        # Get canonical SMILES for each fragment to compare their structure.
        frag_smiles = [Chem.MolToSmiles(frag) for frag in large_frags]
        unique_fragments = set(frag_smiles)
        if len(unique_fragments) >= 2:
            return True, ("Detected a mixture of at least two distinct macromolecular components "
                          "(each with molecular weight >= 100 Da)")
        else:
            return False, ("Multiple fragments were found but they appear to be identical; "
                           "this is more likely a salt than a polymer mixture")
    
    # Note: If the substance is represented as a single fragment,
    # it is usually a defined compound. Even if the compound is very large,
    # by our definition it does not represent a mixture.
    if len(large_frags) == 1:
        frag = large_frags[0]
        if frag.GetNumAtoms() > 100:
            return False, ("A single macromolecule with more than 100 atoms was detected; "
                           "while this may be a large (macromolecular) compound, it is not a mixture "
                           "of different macromolecules as required for a polymer")
    
    # If no large fragments or no mixture of different parts was found, it is not a polymer.
    return False, "Does not meet the criteria for a polymer (mixture of distinct macromolecules)"