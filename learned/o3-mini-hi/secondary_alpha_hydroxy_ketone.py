"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary Alpha-Hydroxy Ketone (Acyloin)

Definition:
    A secondary α‐hydroxy ketone (acyloin) is defined here as a molecule that contains the motif
         R–CH(OH)–C(=O)–R′
    where the “CH(OH)” center is secondary (i.e. it has exactly one hydrogen and three heavy‐atom neighbors)
    and the carbonyl carbon is further attached (via a double bond to oxygen) to at least one carbon (thus not a free acid).

Our strategy is to define a SMARTS query that captures a candidate acyloin motif.
The SMARTS we use is:
    [CH1;X4]([OX2H])-[C;X3](=O)[#6]
This means:
    • [CH1;X4]: an sp3 carbon with exactly one hydrogen
    • ([OX2H]): that carbon is attached to an oxygen that is in an –OH unit
    • -[C;X3](=O): via a single bond the carbon attaches to an sp2 carbon bearing a double bond to an oxygen (a “carbonyl”)
    • [#6]: and the carbonyl carbon is further attached to at least one carbon atom.
    
This query avoids many cases where a free carboxylate or acid appears and (we hope) matches the true acyloin motif.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule contains a secondary alpha-hydroxy ketone (acyloin) motif.
    The motif is defined as R–CH(OH)–C(=O)–R′, where the CH(OH) center is secondary.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the acyloin motif is found, False otherwise.
        str : Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS query.
    # [CH1;X4]: sp3 carbon with exactly one hydrogen.
    # ([OX2H]): attached -OH group.
    # -[C;X3](=O): single bond to an sp2 carbon with a double-bonded oxygen.
    # [#6]: the carbonyl carbon is further bonded to a carbon.
    acyloin_smarts = "[CH1;X4]([OX2H])-[C;X3](=O)[#6]"
    query = Chem.MolFromSmarts(acyloin_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern"
    
    # Check for a substructure match.
    if mol.HasSubstructMatch(query):
        return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    else:
        return False, "Does not contain the required secondary alpha-hydroxy ketone motif"

# Optional: run some tests if executed as a script.
if __name__ == '__main__':
    # Here are a couple tests:
    examples = {
        "(S)-phenylacetylcarbinol": "C=1C=CC=CC1[C@@H](C(C)=O)O",
        "dTDP-4-dehydro-2,6-dideoxy-D-glucose": "C[C@H]1OC(C[C@@H](O)C1=O)OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C[C@@H]1O)n1cc(C)c(=O)[nH]c1=O",
        "2-hydroxyarachidate (false positive expected)": "CCCCCCCCCCCCCCCCCC(O)C([O-])=O",
        "(S)-benzoin": "O[C@H](C(=O)c1ccccc1)c1ccccc1"
    }
    
    for name, smi in examples.items():
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"{name}: {result} -- {reason}")