"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide must contain four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized peptide bond pattern N-C(=O)
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX2H1,NX2H2][C](=O)")  # captures peptide/cyclic bonds
    
    # Checking at least 3 peptide linkages implies a tetrapeptide
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) < 3:
        return False, f"Contains {len(peptide_matches)} peptide bonds, expected at least 3"

    # Refine the amino acid backbone pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX2H1,NX2H2][CX4][CX3](=O)")  # generalized backbone with flexibility
    
    # Look for the presence of four amino residues
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) < 4:
        return False, f"Detected {len(amino_acid_matches)} amino acid residues, expected 4"

    return True, "Contains four amino-acid residues connected by peptide linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48030',
                          'name': 'tetrapeptide',
                          'definition': 'Any molecule that contains four '
                                        'amino-acid residues connected by '
                                        'peptide linkages.',
                          'parents': ['CHEBI:25676'],
                          'xrefs': ['KEGG:C01960'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.42801556420233466 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'N([C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(CF)=O)CC(=O)OC)=O)C(C)C)=O)CCC(OC)=O)=O)CC(OC)=O)C(OCC=1C=CC=CC1)=O '
               'NAME: Z-DEVD-FMK REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CSCC[C@H](NC(=O)[C@H](CCSC)NC(=O)[C@H](CCSC)NC(=O)[C@@H](N)CS)C(O)=O '
               'NAME: Cys-Met-Met-Met REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CO)C(O)=O '
               'NAME: Asp-Leu-Val-Ser REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](N)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O '
               'NAME: Ile-Leu-Val-Tyr REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Asp-Val-Pro-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CN[C@@H](C)C(=O)N[C@@H](C1CCCCC1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](C(c1ccccc1)c1ccccc1)C(=O)OC '
               'NAME: methyl '
               '1-{(2S)-2-cyclohexyl-2-[(N-methyl-L-alanyl)amino]acetyl}-L-prolyl-beta-phenyl-L-phenylalaninate '
               'REASON: CORRECT Contains four amino-acid residues connected by '
               'peptide linkages\n'
               ' * SMILES: '
               'C[C@H](N)C(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Ala-Ala-Pro-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)O)CC(=O)O)CC1=CNC2=C1C=CC=C2)CCSC '
               'NAME: Glu-Met-Trp-Asp REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C[C@@H](O)[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Ala-Thr-Ala-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CCC(O)=O)N[C@H](C(=O)N[C@@H](CC(=O)O)C(N[C@H](C(=O)O)CCCNC(=N)N)=O)CC1=CNC2=C1C=CC=C2 '
               'NAME: Glu-Trp-Asp-Arg REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'N[C@@H](CS)C(=O)N[C@@H](CS)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O '
               'NAME: Cys-Cys-His-His REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](CC(O)=O)C(=O)NCC(O)=O '
               'NAME: Ala-Leu-Asp-Gly REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O '
               'NAME: Ala-Val-Asp-Tyr REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'N[C@@H](CC(O)=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(O)=O '
               'NAME: Asp-Pro-Ser-Ser REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Ala-Leu-Thr-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](C(C)C)C(O)=O '
               'NAME: Asp-Leu-Phe-Val REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](CC(O)=O)C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Ala-Val-Asp-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CS)C(O)=O '
               'NAME: Arg-Phe-Phe-Cys REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CSCC[C@H](NC(=O)[C@@H](N)CS)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O '
               'NAME: Cys-Met-Thr-Tyr REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CC(C)[C@H](NC(=O)[C@@H]1CCCN1)C(=O)NCC(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Pro-Val-Gly-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: C(=O)([C@@H](N)C)NCC(=O)N[C@H](C(=O)NCC(=O)O)C '
               'NAME: Ala-Gly-Ala-Gly REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O '
               'NAME: Ala-Ala-Asp-Asp REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CCCNC(=N)N)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)O)CC(=O)O)[C@H](O)C)[C@H](O)C '
               'NAME: Arg-Thr-Thr-Asp REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'CSCC[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H]([C@@H](C)O)C(=O)N1CCC[C@H]1C(O)=O '
               'NAME: Asp-Met-Thr-Pro REASON: CORRECT Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'C[C@H](N)C(=O)N[C@H](CCC(=O)N[C@@H](CCCCN)C(=O)N[C@H](C)C(O)=O)C(O)=O '
               'NAME: L-alanyl-gamma-D-glutamyl-L-lysyl-D-alanine REASON: '
               'CORRECT Contains four amino-acid residues connected by peptide '
               'linkages\n'
               'False positives: SMILES: '
               'O=C(O)[C@H](NC(=O)C(NC(=O)C(NC(=O)[C@@H](N)C)C(C)C)C(C)C)CC1=CC=C(O)C=C1 '
               'NAME: Aspergillide G REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'OC(=O)CN1C(CN(CCN(CCN(CC1)CC(O)=O)CC(O)=O)CC(O)=O)CCCCNC(=O)CCN2C(=O)C=CC2=O '
               'NAME: '
               '2-[4,7,10-Tris(carboxymethyl)-6-[4-[3-(2,5-dioxopyrrol-1-yl)propanoylamino]butyl]-1,4,7,10-tetrazacyclododec-1-yl]acetic '
               'acid REASON: WRONGLY CLASSIFIED Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'S(CCC(NC(=O)C(NC(=O)C(N)C)CCSC)C(=O)NC(CC1=CC=C(O)C=C1)C(O)=O)C '
               'NAME: Ala-Met-Met-Tyr REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(NC(C(C)C)C(=O)NC(CO)C(O)=O)C(NC(=O)C(N)CCCN=C(N)N)CCCCN '
               'NAME: Arg-Lys-Val-Ser REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'C[C@@H](O)[C@H](NC(=O)[C@H](CCCNC(N)=N)NC(=O)c1cccc(O)c1O)C(=O)O[C@H](C)[C@H](NC(=O)[C@H](CCCNC(N)=N)NC(=O)c1cccc(O)c1O)C(=O)O[C@H](C)[C@H](NC(=O)[C@H](CCCNC(N)=N)NC(=O)c1cccc(O)c1O)C(O)=O '
               'NAME: Tribenarthin REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'C[C@H]1OC(=O)[C@@H](NC(=O)CNC(=O)c2cccc(O)c2O)[C@@H](C)OC(=O)[C@@H](NC(=O)CNC(=O)c2cccc(O)c2O)[C@@H](C)OC(=O)[C@H]1NC(=O)CNC(=O)c1cccc(O)c1O '
               'NAME: corynebactin REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N1C(CCC1)C(=O)NC(CC=2NC=NC2)C(O)=O)C(NC(=O)C(N)CC(O)=O)CC3=CC=CC=C3 '
               'NAME: Asp-Phe-Pro-His REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](N)CC1=CC=CC=C1)C(C)C)C(C)C)CC2=CC=C(O)C=C2 '
               'NAME: D-Phe-L-Val-D-Val-L-Tyr REASON: WRONGLY CLASSIFIED '
               'Contains four amino-acid residues connected by peptide '
               'linkages\n'
               ' * SMILES: '
               'O=C1N(C=CC(N1)=O)C2O/C(=C\\NC(=O)C(NC(=O)C(NC(=O)NC(C(=O)O)CC=3C4=C(C=CC=C4)NC3)C)C(N(C(=O)C(N)C)C)C)/CC2O '
               'NAME: Pacidamycin D REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C1N(C=CC(N1)=O)[C@@H]2O/C(=C\\NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)N[C@H](C(=O)O)CC3=CC=CC=C3)C)[C@@H](N(C(=O)[C@@H](N)CC4=CC(O)=CC=C4)C)C)/C[C@H]2O '
               'NAME: Pacidamycin 5 REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'OC(C(NC(=O)CNC(=O)CN)C(=O)NC(CC1=CC=C(O)C=C1)C(O)=O)C NAME: '
               'Gly Gly Thr Tyr REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(NC(CC1=CC=CC=C1)C(=O)NC(CC2=CC=CC=C2)C(=O)NC(CC=3C=4C(NC3)=CC=CC4)C(O)=O)C(N)C(CC)C '
               'NAME: Ile-Phe-Phe-Trp REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'S(C[C@H](N(C(=O)[C@@H](N(C(=O)CNC(=O)[C@H](NC(=O)C1=NC=2C(=CC=CC2)C=C1O)CSC)C)CSC)C)C(=O)OC)C '
               'NAME: Thiochondrilline A REASON: WRONGLY CLASSIFIED Contains '
               'four amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N1[C@@H](CCC1)C(=O)N[C@@H](CC2=CC=C(O)C=C2)C(=O)N[C@H]([C@@H](CC)C)C(O)=O)[C@H](N)[C@@H](CC)C '
               'NAME: Neuromedin N (1-4) REASON: WRONGLY CLASSIFIED Contains '
               'four amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N1C(CCC1)C(O)=O)C(NC(=O)C(NC(=O)C(N)C(CC)C)C(CC)C)C(CC)C '
               'NAME: Ile-Ile-Ile-Pro REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C1N[C@@H]2C(=O)N[C@H](C(=O)O)CC3=CC=C(OC=4C=C(C2)C=C(C5=C(OC)C=CC(C[C@@H](C(N[C@H]1[C@H](CC)C)=O)N)=C5)C4O)C=C3 '
               'NAME: Cittilin A REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'ClCCCCCCCC(N)CNC(C(=O)N1C(C(=O)NC(C(=O)NC(C(=O)O)CC2=CC=C(O)C=C2)CC3=CC=C(O)C=C3)CCC1)C(C)C '
               'NAME: Microginin 715 REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C1NC(C(=O)NCCCCC(C(NC1C(CC)C)=O)NC(=O)NC(C(=O)O)C(CC)C)CCC2=CC=CC=C2 '
               'NAME: Namalide D REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: O=C(NCC(=O)NC(C)C(O)=O)C(NC(=O)C(N)C)C NAME: Ala '
               'Ala Gly Ala REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'CCC(C)C(C(=O)NC(CCCN=C(N)N)C(=O)NC(CC(C)C)C(=O)O)NC(=O)C(CC1=CN=CN1)N '
               'NAME: '
               '2-[[2-[[2-[[2-amino-3-(1H-imidazol-5-yl)-1-oxopropyl]amino]-3-methyl-1-oxopentyl]amino]-5-(diaminomethylideneamino)-1-oxopentyl]amino]-4-methylpentanoic '
               'acid REASON: WRONGLY CLASSIFIED Contains four amino-acid '
               'residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N[C@H](C(C)C)C(=O)N[C@@H]([C@@H](O)C)C(O)=O)[C@@H](NC(=O)[C@H](N)CCCCN)C(C)C '
               'NAME: Dynorphin B (10-13) REASON: WRONGLY CLASSIFIED Contains '
               'four amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)O)[C@@H](CC)C)CC(C)C)C)C)CCO)C1=CC=CC=C1 '
               'NAME: Trichopeptide A REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: SCC(NC(=O)C(NC(=O)C(N)C)C(O)C)C(=O)NC(CS)C(O)=O '
               'NAME: Ala-Thr-Cys-Cys REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C1OC[C@@H]([C@H](O)[C@@H](C(=O)N[C@@H](C(N[C@H](C(N[C@H]1CO)=O)C(=O)C=2C3=C(C=CC=C3)NC2)=O)CCC(=O)O)C)CCCCCCCC '
               'NAME: 15G256delta REASON: WRONGLY CLASSIFIED Contains four '
               'amino-acid residues connected by peptide linkages\n'
               ' * SMILES: '
               'O=C(N[C@H](CCC(=O)N[C@H](CCC=O)C(O)=O)C(O)=O)CC[C@@H](NC(=O)CC[C@@H](N)C(O)=O)C(O)=O '
               'NAME: Poly-g-D-glutamate REASON: WRONGLY CLASSIFIED Contains '
               'four amino-acid residues connected by peptide linkages\n'
               'False negatives: SMILES: '
               'N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(O)=O '
               'NAME: Glu-Asp-Gln-Gln REASON: MISSED Contains 5 peptide bonds, '
               'expected 3\n'
               ' * SMILES: '
               'CSCC[C@H](NC(=O)[C@H]1N(Cc2ccccc12)C(=O)[C@@H](NC(=O)[C@H](CC1CCCCC1)NC(C)=O)C1CCN(CC1)C(N)=N)C(N)=O '
               'NAME: Ac-Cha-Gpg-Disc-Met-NH2 REASON: MISSED Contains 5 '
               'peptide bonds, expected 3\n'
               ' * SMILES: '
               'CC(C)C[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCC(N)=O)C(O)=O '
               'NAME: Ala-Leu-Thr-Gln REASON: MISSED Contains 4 peptide bonds, '
               'expected 3\n'
               ' * SMILES: '
               'NCCCC[C@H](NC(=O)[C@@H](N)CC(N)=O)C(=O)N[C@@H](CC1=CN=CN1)C(=O)N[C@@H](CC1=CN=CN1)C(O)=O '
               'NAME: Asn-Lys-His-His REASON: MISSED Contains 4 peptide bonds, '
               'expected 3\n'
               ' * SMILES: '
               'NC(CCCCNC(=O)CCCCCNC(=O)c1ccc(o1)[N+]([O-])=O)C(=O)NCCCCC(N)C(=O)NCCCCC(N)C(=O)NC(CCCCNC(=O)CCCc1c[nH]c2ccccc12)C(N)=O '
               'NAME: '
               'N(6)-[N(6)-(N(6)-{6-[(5-nitro-2-furoyl)amino]hexanoyl}lysyl)lysyl]lysyl-N(6)-[4-(indol-3-yl)butanoyl]lysinamide '
               'REASON: MISSED Contains 6 peptide bonds, expected 3\n'
               ' * SMILES: '
               'CCCC[C@H](NC(=O)[C@@H]1Cc2ccccc2CN1C(=O)[C@@H](NC(=O)[C@H](CC1CCCCC1)NC(C)=O)C1CCN(CC1)C(N)=N)C(=O)N(C)Cc1ccccc1 '
               'NAME: Ac-Cha-Gpg-Tic-Nle-NMeBn REASON: MISSED Contains 5 '
               'peptide bonds, expected 3\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@H](OC[C@H](NC(=O)CNC(C)=O)C(=O)NCC(=O)N[C@@H](OC(C)C)C(O)=O)[C@@H](O)[C@@H](O)[C@@H]1O)[C@H](O)CO '
               'NAME: AcGly-Ser(O-L-alpha-D-Hep)-Gly-Leu REASON: MISSED '
               'Contains 4 peptide bonds, expected 3\n'
               ' * SMILES: '
               'N[C@@H](CC(O)=O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CCC(N)=O)C(O)=O '
               'NAME: Asp-Phe-Asp-Gln REASON: MISSED Contains 4 peptide bonds, '
               'expected 3\n'
               ' * SMILES: '
               'N[C@@H](CCCNC(N)=N)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCCNC(N)=N)C(O)=O '
               'NAME: Arg-Asn-Gln-Arg REASON: MISSED Contains 5 peptide bonds, '
               'expected 3\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@H](OC[C@H](NC(=O)CNC(C)=O)C(=O)NCC(=O)N[C@@H](OC(C)C)C(O)=O)[C@@H](O)[C@@H](O)[C@@H]1O)[C@@H](O)CO '
               'NAME: AcGly-Ser(O-D-alpha-D-Hep)-Gly-Leu REASON: MISSED '
               'Contains 4 peptide bonds, expected 3\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'Contains 2 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'Contains 0 peptide bonds, '
                                               'expected at least 3'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'Contains 1 peptide bonds, '
                                               'expected at least 3'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]([C@@H](C)CC)([C@@H](CC(=O)N1CCC[C@@]1([H])[C@H](OC)[C@@H](C)C(=O)NCCC1=CC=CC=C1)OC)N(C)C(=O)[C@@H](NC(=O)[C@H](C(C)C)N(C)C)C(C)C',
                                      'name': 'soblidotin',
                                      'reason': 'Detected 2 amino acid '
                                                'residues, expected 4'},
                                  {   'smiles': 'CC(C)C[C@@H](N)[C@H](O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(O)=O',
                                      'name': 'amastatin',
                                      'reason': 'Detected 3 amino acid '
                                                'residues, expected 4'},
                                  {   'smiles': 'CC(=O)N[C@H](Cc1ccccc1)C(=O)NCCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N1CCC[C@@H]1C(N)=O',
                                      'name': 'Ac-D-Phe-betaAla-L-His-D-Pro-NH2',
                                      'reason': 'Detected 3 amino acid '
                                                'residues, expected 4'},
                                  {   'smiles': 'NC(=O)CCOCCOCCOCCOCCOCCNC(=O)C(CCCCNC(=O)c1ccc(o1)[N+]([O-])=O)NC(=O)CCOCCOCCNC(=O)C1CCN(CC1)C(=O)C1CCN(CC1)C(=O)C1CCN(CC1)C(=O)C1CCN(CC1)C(=O)CCOCCOCCNC(=O)CCCc1c[nH]c2ccccc12',
                                      'name': 'N-(18-amino-18-oxo-3,6,9,12,15-pentaoxaoctadec-1-yl)-N(2)-{3-[2-(2-{[(1-{[1-({1-[(1-{3-[2-(2-{[4-(indol-3-yl)butanoyl]amino}ethoxy)ethoxy]propanoyl}piperidin-4-yl)carbonyl]piperidin-4-yl}carbonyl)piperidin-4-yl]carbonyl}piperidin-4-yl)carbonyl]amino}ethoxy)ethoxy]propanoyl}-N(6)-(5-nitro-2-furoyl)lysinamide',
                                      'reason': 'Detected 1 amino acid '
                                                'residues, expected 4'},
                                  {   'smiles': 'O=C(N[C@H](C(N[C@H](C(NCCC(=O)O)=O)CCCCNC(CCCCC1SCC2C1NC(N2)=O)=O)=O)CCCCN)[C@@H](N(C3=CC=C(C=4C3=NON4)[N+]([O-])=O)C)CC(=O)O',
                                      'name': 'NBD-N(Me)-DKK(Biotin)betaA-OH',
                                      'reason': 'Detected 3 amino acid '
                                                'residues, expected 4'}],
    'attempt': 4,
    'reasoning': 'To improve the classification of tetrapeptides based on '
                 'their SMILES strings, we need to address some issues '
                 'observed in the previous attempt:\n'
                 '\n'
                 '1. **Pattern Definition Errors**: The peptide bond pattern '
                 '`[N][C;$(C(=O)[N])]C` might not be accurately capturing the '
                 'variety of peptide bonds found in the examples. We need a '
                 'more generalized pattern to correctly identify peptide '
                 'linkages.\n'
                 '\n'
                 '2. **Counting Peptide Bonds**: The condition '
                 '`len(peptide_matches) != 3` for a tetrapeptide may not be '
                 'flexible enough considering potential variations in SMILES '
                 'representations, such as internal cyclic structures or '
                 'protecting groups which mask true peptide bonds in a simple '
                 'count.\n'
                 '\n'
                 '3. **Amino Acid Residue Representation**: The SMARTS pattern '
                 'to identify amino acid residues might need refinement. '
                 'Besides, a tetrapeptide specifically requires four residues, '
                 'but our pattern might not be robust enough due to '
                 'stereochemistry and potential protecting groups.\n'
                 '\n'
                 '**Hypothesis for Improvement**:\n'
                 '- Create a flexible pattern to recognize a variety of '
                 'peptide linkages.\n'
                 '- Improve amino acid backbone pattern to account for chiral '
                 'centers and potential extensions such as side-chain '
                 'protections or modifications.\n'
                 '- Ensure the conditions are flexible enough to accommodate '
                 'different SMILES representations while accurately '
                 'identifying tetrapeptides.\n'
                 '\n'
                 'Here is a revised implementation with these considerations:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 97,
    'num_false_positives': 100,
    'num_true_negatives': 4304,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.49238578680203043,
    'recall': 0.9509803921568627,
    'f1': 0.6488294314381271,
    'accuracy': 0.9766977363515313,
    'negative_predictive_value': 0.9988396379670457}