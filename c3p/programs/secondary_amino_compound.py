"""
Classifies: CHEBI:50995 secondary amino compound
"""
from rdkit import Chem

def is_secondary_amino_compound(smiles: str):
    """
    Determines if a molecule is a secondary amino compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amino compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms in the molecule
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    
    for nitrogen in nitrogen_atoms:
        # Check if the nitrogen has exactly two carbon neighbors and no hydrogen
        carbon_neighbors = [neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetSymbol() == 'C']
        hydrogen_neighbors = [neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetSymbol() == 'H']
        
        if len(carbon_neighbors) == 2 and len(hydrogen_neighbors) == 0:
            return True, "Secondary amino compound found"
    
    return False, "No secondary amino compound found"

# Test cases
smiles_list = [
    "[H][C@]12CN(C)[C@]([H])(CN1)CC1=C[C@@]([H])(C(=O)CC1)[C@@]1([H])C=C(CCC1=O)C2",
    "C1(=C(CC2=CN=CC=C2)C3=C(C(C)=C1O)SC(NC)=N3)C",
    "O(C1=C2C(=C(C=C1OC)NC(CCCN)C)N=C(C=C2C)OC)C3=CC(C(F)(F)F)=CC=C3",
    "CC[C@H](Nc1ncnc(C)c1Cl)c1ccc(OC(F)F)cc1",
    "Clc1ccccc1Nc1nc(Cl)nc(Cl)n1",
    "OCc1cc(ccc1O)[C@@H](O)CNCCCCCCOCCCCc1ccccc1",
    "Cc1ccc(CNc2ccccc2C(O)=O)cc1",
    "COC1=CC=CC(=C1)NC2=NC=NC3=CC=CC=C32",
    "C1(=C(Cl)C=CC=C1)[C@@H](CNC(C)C)O",
    "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)O)N)=O",
    "[H][C@@]12CCC=C(C)[C@@]1(C)CC[C@](C)(CC1=C(O)C(=O)C=C(NCCS(O)(=O)=O)C1=O)[C@@H]2C",
    "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1",
    "N(C1CCCC1)CC",
    "OC(=O)CCCCCCNC1c2ccccc2CCc2ccccc12",
    "CC(C)NCC(O)COc1cccc2[nH]ccc12",
    "C=1C(=CC=CC1Cl)[C@H](CNCCNC=2C=C(C=CC2)C=3C=CC=C(C3)C(=O)O)O",
    "CC(C)N[C@@H](C)[C@@H](O)COc1ccc(C)c2CCCc12",
    "COc1ccc(Cl)c(Nc2ncnc3cc(OCCCN4CCOCC4)c(OC)cc23)c1",
    "C[C@@H]1CN(CCN1)c1cc2n(cc(C(O)=O)c(=O)c2cc1F)-c1ccc(F)cc1F",
    "OS(=O)(=O)c1cccc(c1)N=Nc1ccc(Nc2ccccc2)cc1",
    "[C@@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1CO)O)O)O)N[C@@H]2[C@@H]([C@H]([C@@H](C(=C2)CO)O)O)O)O",
    "O(CCCCCCCCCC)CCNCCOCCCCCCCCCC",
    "CN[C@@H]1Cc2c[nH]c3cccc([C@H]1\C=C(/C)C=O)c23",
    "CC(C)(C)NC[C@@H](O)COc1cccc2C[C@@H](O)[C@@H](O)Cc12",
    "CC(CCCC(NC(C)C)C)C.Cl",
    "Cc1cn(cn1)-c1cc(NC(=O)c2ccc(C)c(Nc3nccc(n3)-c3cccnc3)c2)cc(c1)C(F)(F)F",
    "BrC=1C=C(F)C(NC=2C(F)=C3N=CN(C3=CC2C(=O)NOCCO)C)=CC1",
    "CNCCc1ccccn1",
    "C1CNCC(C2=C1C=C(C(=C2)O)O)C3=CC=CC=C3",
    "Cc1ccc(cc1)C(=O)Oc1ccc(cc1OC(=O)c1ccc(C)cc1)[C@H](O)CNC(C)(C)C",
    "CC(CCC1=CC=CC=C1)NCC(O)C1=CC(C(N)=O)=C(O)C=C1",
    "OC1=CC=C(\C=C\C(=O)NCCCNCCCCNC(=O)\C=C\C2=CC(O)=C(O)C=C2)C=C1O",
    "NC(=O)CC[C@@H]1NC[C@@]2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O",
    "C[C@@H]1COC(NC2=CC3=C(C=C2)N=CN=C3NC2=CC(Cl)=C(OCC3=NC=CS3)C=C2)=N1",
    "CCNC(CC(N)=O)C(O)=O",
    "C=1C(=C(C2=C(C1NC(=O)C(C)=CC=C[C@@H]([C@@H](OC(N)=O)C(C)=C[C@@H]([C@H]([C@H](C[C@@H](C2)C)OC)O)C)OC)O)NCC=C)O",
    "C=1(C(NC=2C=C(C=CC2)CO)=CC=NC1)S(NC(NC(C)C)=O)(=O)=O",
    "C=12C(=C(N=CN1)NC3CN(CC3)CC=4C=CC=CC4)C(=C(S2)C)C",
    "CC[C@H]1N(C(C)C)C2=NC(NC3=C(OC)C=C(C=C3)C(=O)N[C@H]3CC[C@@H](CC3)N3CCN(CC4CC4)CC3)=NC=C2N(C)C1=O",
    "S(CCNC1CCCCC1)(O)(=O)=O",
    "CC1=CC=CC(C2=NC3=CC=CC=C3C(NC4=CC=NC=C4)=N2)=N1",
    "C(CNC(C1=CC=CC=C1)C1=CC=CC=C1)NC(C1=CC=CC=C1)C1=CC=CC=C1",
    "Cc1ccc(cc1)C(=O)Oc1ccc(cc1OC(=O)c1ccc(C)cc1)C(O)CNC(C)(C)C",
    "C1N(C1)CCNCCN",
    "CNC1=CC(OC)=C(C=C1Cl)C(=O)N[C@@H]1CCN(CC2=CC=CC=C2)[C@@H]1C",
    "C(CCCC(C)C)(NCCC(C)C)C",
    "COC1=C(NC2=NC=C(Cl)C(=N2)C2=CNC3=CC=CC=C23)C=CC(=C1)N1CCC(N)CC1",
    "N(CC(C)C)CCC",
    "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1",
    "C1=C(C=CC(=C1)[C@H](C(CC)CC)N2C=NC=N2)NC=3SC4=C(N3)C=CC=C4",
    "CNC1=CC(=O)C2=C(C=NC(C)=C2C)C1=O",
    "C=12N=C(N=CC1C3=C(N2[C@@H]4CC[C@H](CC4)C)C=NC=C3)NC=5C=CC6=C(N5)CCN(C6)C(CO)=O",
    "C[C@@H](NCCCc1cccc(c1)C(F)(F)F)c1cccc2ccccc12",
    "CN1C=C(CNCC2CCN(CC2)C2=NC=C(C=N2)C(=O)NO)C2=CC=CC=C12",
    "N(CCCCNCCCN)CCCNCCCCN",
    "CNC[C@H](O)c1cccc(O)c1",
    "CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(c1)C(C)=O"
]

for smiles in smiles_list:
    result, reason = is_secondary_amino_compound(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50995',
                          'name': 'secondary amino compound',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing two hydrogen '
                                        'atoms by organyl groups.',
                          'parents': ['CHEBI:50047']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              '[H][C@]12CN(C)[C@]([H])(CN1)CC1=C[C@@]([H])(C(=O)CC1)[C@@]1([H])C=C(CCC1=O)C2\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C1(=C(CC2=CN=CC=C2)C3=C(C(C)=C1O)SC(NC)=N3)C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'O(C1=C2C(=C(C=C1OC)NC(CCCN)C)N=C(C=C2C)OC)C3=CC(C(F)(F)F)=CC=C3\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC[C@H](Nc1ncnc(C)c1Cl)c1ccc(OC(F)F)cc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: Clc1ccccc1Nc1nc(Cl)nc(Cl)n1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: OCc1cc(ccc1O)[C@@H](O)CNCCCCCCOCCCCc1ccccc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: Cc1ccc(CNc2ccccc2C(O)=O)cc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: COC1=CC=CC(=C1)NC2=NC=NC3=CC=CC=C32\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C1(=C(Cl)C=CC=C1)[C@@H](CNC(C)C)O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              '[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)O)N)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              '[H][C@@]12CCC=C(C)[C@@]1(C)CC[C@](C)(CC1=C(O)C(=O)C=C(NCCS(O)(=O)=O)C1=O)[C@@H]2C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: N(C1CCCC1)CC\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: OC(=O)CCCCCCNC1c2ccccc2CCc2ccccc12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(C)NCC(O)COc1cccc2[nH]ccc12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C=1C(=CC=CC1Cl)[C@H](CNCCNC=2C=C(C=CC2)C=3C=CC=C(C3)C(=O)O)O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(C)N[C@@H](C)[C@@H](O)COc1ccc(C)c2CCCc12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: COc1ccc(Cl)c(Nc2ncnc3cc(OCCCN4CCOCC4)c(OC)cc23)c1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C[C@@H]1CN(CCN1)c1cc2n(cc(C(O)=O)c(=O)c2cc1F)-c1ccc(F)cc1F\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: OS(=O)(=O)c1cccc(c1)N=Nc1ccc(Nc2ccccc2)cc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              '[C@@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1CO)O)O)O)N[C@@H]2[C@@H]([C@H]([C@@H](C(=C2)CO)O)O)O)O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: O(CCCCCCCCCC)CCNCCOCCCCCCCCCC\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CN[C@@H]1Cc2c[nH]c3cccc([C@H]1\\C=C(/C)C=O)c23\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(C)(C)NC[C@@H](O)COc1cccc2C[C@@H](O)[C@@H](O)Cc12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(CCCC(NC(C)C)C)C.Cl\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'Cc1cn(cn1)-c1cc(NC(=O)c2ccc(C)c(Nc3nccc(n3)-c3cccnc3)c2)cc(c1)C(F)(F)F\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: BrC=1C=C(F)C(NC=2C(F)=C3N=CN(C3=CC2C(=O)NOCCO)C)=CC1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CNCCc1ccccn1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C1CNCC(C2=C1C=C(C(=C2)O)O)C3=CC=CC=C3\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'Cc1ccc(cc1)C(=O)Oc1ccc(cc1OC(=O)c1ccc(C)cc1)[C@H](O)CNC(C)(C)C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(CCC1=CC=CC=C1)NCC(O)C1=CC(C(N)=O)=C(O)C=C1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'OC1=CC=C(\\C=C\\C(=O)NCCCNCCCCNC(=O)\\C=C\\C2=CC(O)=C(O)C=C2)C=C1O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'NC(=O)CC[C@@H]1NC[C@@]2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C[C@@H]1COC(NC2=CC3=C(C=C2)N=CN=C3NC2=CC(Cl)=C(OCC3=NC=CS3)C=C2)=N1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CCNC(CC(N)=O)C(O)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C=1C(=C(C2=C(C1NC(=O)C(C)=CC=C[C@@H]([C@@H](OC(N)=O)C(C)=C[C@@H]([C@H]([C@H](C[C@@H](C2)C)OC)O)C)OC)O)NCC=C)O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C=1(C(NC=2C=C(C=CC2)CO)=CC=NC1)S(NC(NC(C)C)=O)(=O)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C=12C(=C(N=CN1)NC3CN(CC3)CC=4C=CC=CC4)C(=C(S2)C)C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'CC[C@H]1N(C(C)C)C2=NC(NC3=C(OC)C=C(C=C3)C(=O)N[C@H]3CC[C@@H](CC3)N3CCN(CC4CC4)CC3)=NC=C2N(C)C1=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: S(CCNC1CCCCC1)(O)(=O)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC1=CC=CC(C2=NC3=CC=CC=C3C(NC4=CC=NC=C4)=N2)=N1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C(CNC(C1=CC=CC=C1)C1=CC=CC=C1)NC(C1=CC=CC=C1)C1=CC=CC=C1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'Cc1ccc(cc1)C(=O)Oc1ccc(cc1OC(=O)c1ccc(C)cc1)C(O)CNC(C)(C)C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C1N(C1)CCNCCN\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'CNC1=CC(OC)=C(C=C1Cl)C(=O)N[C@@H]1CCN(CC2=CC=CC=C2)[C@@H]1C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C(CCCC(C)C)(NCCC(C)C)C\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'COC1=C(NC2=NC=C(Cl)C(=N2)C2=CNC3=CC=CC=C23)C=CC(=C1)N1CCC(N)CC1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: N(CC(C)C)CCC\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C1=C(C=CC(=C1)[C@H](C(CC)CC)N2C=NC=N2)NC=3SC4=C(N3)C=CC=C4\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CNC1=CC(=O)C2=C(C=NC(C)=C2C)C1=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: '
              'C=12N=C(N=CC1C3=C(N2[C@@H]4CC[C@H](CC4)C)C=NC=C3)NC=5C=CC6=C(N5)CCN(C6)C(CO)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: C[C@@H](NCCCc1cccc(c1)C(F)(F)F)c1cccc2ccccc12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CN1C=C(CNCC2CCN(CC2)C2=NC=C(C=N2)C(=O)NO)C2=CC=CC=C12\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: N(CCCCNCCCN)CCCNCCCCN\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CNC[C@H](O)c1cccc(O)c1\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n'
              'SMILES: CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(c1)C(C)=O\n'
              'Result: True\n'
              'Reason: Secondary amino compound found\n'
              '\n',
    'num_true_positives': 57,
    'num_false_positives': 7,
    'num_true_negatives': 13,
    'num_false_negatives': 0,
    'precision': 0.890625,
    'recall': 1.0,
    'f1': 0.9421487603305785,
    'accuracy': None}