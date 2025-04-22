# Curated Validation

In order to explore further the extent to which incorrect benchmarks may be detrimental to performance and under-estimating true performance, we systematically investigated a set of individual structure-class pairs. We took  all chemical classes from the ensemble model with a program that scored 0.8 or higher. For each of these classes, we sampled 3 false positive and 3 false negative structures. We ran these through ChEBACR (see methods) to get a second opinion on 221 cases in which both C3PO and multiple LLMs were in agreement.

* [this sheet](https://docs.google.com/spreadsheets/d/1EifaWOJb-ClvPGK4GW7dIOhD4WB9XLAtdaBpDGqjvGI/edit?gid=298616687#gid=298616687)
* [as xlsx](chebacr-chebi-curators.xslx)

## Example: octadecanedioic acid

class ID | CHEBI:18310
---|---
class name | alkane
structure  | C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O
structure name | octadecanedioic acid
gpt-4o response (with literature search via citeseek) | NO, Alkanes are hydrocarbons, meaning they consist solely of hydrogen and carbon atoms, with the general formula CnH2n+2. They are saturated, containing only single bonds between carbon atoms. Octadecanedioic acid, on the other hand, contains not only carbon and hydrogen but also oxygen atoms due to the presence of carboxylic acid groups (–COOH), which means it is not a pure hydrocarbon. The structure given in the SMILES representation indicates that octadecanedioic acid is a dicarboxylic acid, with two carboxylic acid functional groups. Therefore, it cannot be classified as an alkane.
Claude response | "NO. Alkanes are hydrocarbons that contain only single bonds between carbon atoms and are fully saturated with hydrogen atoms. The general formula for alkanes is CnH2n+2. However, octadecanedioic acid is not an alkane because it contains two carboxyl groups (-COOH) at the opposite ends of the carbon chain. The presence of these carboxyl groups, which contain double-bonded oxygen atoms, means that the compound is not fully saturated with hydrogen atoms and is not an alkane. Instead, octadecanedioic acid is classified as a dicarboxylic acid, which is a type of carboxylic acid containing two carboxyl groups. It belongs to the broader class of organic compounds known as oxoacids, which contain one or more oxygen atoms and at least one hydroxyl group (-OH) attached to a carbon atom. Therefore, the correct classification for octadecanedioic acid is a dicarboxylic acid, not an alkane."
DeepSeek chat response | NO. Octadecanedioic acid has the SMILES 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O', which indicates the presence of two carboxyl groups (-COOH) at the ends of the molecule. Alkanes, by definition, are saturated hydrocarbons consisting solely of hydrogen and carbon atoms, with no functional groups. The presence of carboxyl groups in octadecanedioic acid disqualifies it from being classified as an alkane. Instead, it is a dicarboxylic acid.
gpt-3.5 response | "NO Explanation: Octadecanedioic acid, with the chemical formula C18H34O4, contains oxygen atoms in addition to carbon and hydrogen atoms, which makes it a carboxylic acid and not an alkane. Alkanes do not contain oxygen atoms and have the general formula CnH2n+2."
gpt-4o response | NO, the statement is not true. Octadecanedioic acid, with the SMILES notation 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O', is not an alkane. Alkanes are acyclic hydrocarbons with the general formula CnH2n+2 and only consist of carbon and hydrogen atoms, with all carbon-carbon bonds being single bonds. Octadecanedioic acid, on the other hand, contains two carboxylic acid groups (indicated by the '(=O)O' in the SMILES notation) which means it has carbon-oxygen double bonds and hydroxyl groups. These functional groups make it a dicarboxylic acid rather than a purely hydrocarbon alkane. Therefore, it does not fit the definition of an alkane.
Curator: response | NO, is a dicarboxylic acid and a fatty acid
C3P: fatty acid | Yes (0.7) found free acid group (\u2013C(=O)OH) with a contiguous aliphatic chain of 18 carbons starting at the carbonyl carbon
C3P: dicarboxylic acid | Yes (0.25) Molecule contains exactly two free carboxyl groups
C3P: alkane | No (0.99) Molecule contains non-carbon and non-hydrogen atoms
