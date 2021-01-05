let HBB_GENE =
  "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGA" +
  "GAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTG" +
  "GTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGG" +
  "TTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTT" +
  "GGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAA" +
  "GGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTT" +
  "GCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACG" +
  "CTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGGATAAGTAACAGGGTACA" +
  "GTTTAGAATGGGAAACAGACGAATGATTGCATCAGTGTGGAAGTCTCAGGATCGTTTTAGTTTCTTTTATTTG" +
  "CTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATT" +
  "ATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAA" +
  "CTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTA" +
  "CTTTATTTTCTTTTATTTTTAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAAT" +
  "ATGTGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATA" +
  "TACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTAT" +
  "CATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATCTCTGCAT" +
  "ATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCT" +
  "ACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAAT" +
  "CATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTT" +
  "GGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCC" +
  "ACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTA" +
  "CTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC";

const START_CODON = "AUG";
const STOP_CODON = "STOP";

const sickleCellAnemia = (gene) => {
  gene = gene.split("");
  gene.splice(69, 1, "U");
  return gene.join("");
};

const intronScanner = {
  start: /GGU(?:G|U)(?:G|A)GU/g,
  end: /CCC(?:U|A)(?:U|C)AG/g,
};

const AMINO_ACIDS = {
  W: { key: "-trypto-*", name: "Tryptophan" },
  C: { key: "-cyst-*", name: "Cystine" },
  G: { key: "-gly-*", name: "Glycine" },
  R: { key: "-arg-*", name: "Arginine" },
  S: { key: "-ser-*", name: "Serine" },
  T: { key: "-threo-*", name: "Threonine" },
  A: { key: "-ala-*", name: "Alanine" },
  P: { key: "-prol-*", name: "Proline" },
  F: { key: "-phenyl-*", name: "Phenylalanine" },
  L: { key: "-leu-*", name: "Leucine" },
  V: { key: "-val-*", name: "Valine" },
  I: { key: "-iso-*", name: "Isoleucine" },
  M: { key: "-meth-*", name: "Methionine" },
  Q: { key: "-glut-*", name: "Glutamine" },
  H: { key: "-hist-*", name: "Histidine" },
  D: { key: "-asp-*", name: "Aspartic acid" },
  E: { key: "-glutacid-*", name: "Glutamic acid" },
  K: { key: "-lys-*", name: "Lysine" },
  N: { key: "-aspa-*", name: "Asparagine" },
  Y: { key: "-tyro-*", name: "Tyrosine" },
};

const CODONS = {
  UUU: AMINO_ACIDS.F,
  UUC: AMINO_ACIDS.F,
  UUA: AMINO_ACIDS.L,
  UUG: AMINO_ACIDS.L,
  UCU: AMINO_ACIDS.S,
  UCC: AMINO_ACIDS.S,
  UCA: AMINO_ACIDS.S,
  UCG: AMINO_ACIDS.S,
  UAU: AMINO_ACIDS.Y,
  UAC: AMINO_ACIDS.Y,
  UAA: STOP_CODON,
  UAG: STOP_CODON,
  UGU: AMINO_ACIDS.C,
  UGC: AMINO_ACIDS.C,
  UGA: STOP_CODON,
  UGG: AMINO_ACIDS.W,
  CUU: AMINO_ACIDS.L,
  CUC: AMINO_ACIDS.L,
  CUA: AMINO_ACIDS.L,
  CUG: AMINO_ACIDS.L,
  CCU: AMINO_ACIDS.P,
  CCC: AMINO_ACIDS.P,
  CCA: AMINO_ACIDS.P,
  CCG: AMINO_ACIDS.P,
  CAU: AMINO_ACIDS.H,
  CAC: AMINO_ACIDS.H,
  CAA: AMINO_ACIDS.Q,
  CAG: AMINO_ACIDS.Q,
  CGU: AMINO_ACIDS.R,
  CGC: AMINO_ACIDS.R,
  CGA: AMINO_ACIDS.R,
  CGG: AMINO_ACIDS.R,
  AUU: AMINO_ACIDS.I,
  AUC: AMINO_ACIDS.I,
  AUA: AMINO_ACIDS.I,
  AUG: AMINO_ACIDS.M,
  ACU: AMINO_ACIDS.T,
  ACC: AMINO_ACIDS.T,
  ACA: AMINO_ACIDS.T,
  ACG: AMINO_ACIDS.T,
  AAU: AMINO_ACIDS.N,
  AAC: AMINO_ACIDS.N,
  AAA: AMINO_ACIDS.K,
  AAG: AMINO_ACIDS.K,
  AGU: AMINO_ACIDS.S,
  AGC: AMINO_ACIDS.S,
  AGA: AMINO_ACIDS.R,
  AGG: AMINO_ACIDS.R,
  GUU: AMINO_ACIDS.V,
  GUC: AMINO_ACIDS.V,
  GUA: AMINO_ACIDS.V,
  GUG: AMINO_ACIDS.V,
  GCU: AMINO_ACIDS.A,
  GCC: AMINO_ACIDS.A,
  GCA: AMINO_ACIDS.A,
  GCG: AMINO_ACIDS.A,
  GAU: AMINO_ACIDS.D,
  GAC: AMINO_ACIDS.D,
  GAA: AMINO_ACIDS.E,
  GAG: AMINO_ACIDS.E,
  GGU: AMINO_ACIDS.G,
  GGC: AMINO_ACIDS.G,
  GGA: AMINO_ACIDS.G,
  GGG: AMINO_ACIDS.G,
};

const rnaPolymerase = (nucleotide) => {
  return nucleotide.replace(/T/g, "U");
};

const splicing = (rnaConverted) => {
  let start = rnaConverted.match(intronScanner.start);
  let end = rnaConverted.match(intronScanner.end);

  if (!start || !end || start.length !== end.length) {
    console.error("Splicing Error");
    return rnaConverted;
  }

  let counter = 0;
  do {
    let intronStart = rnaConverted.indexOf(start[counter]) + 1;
    let intronEnd = rnaConverted.indexOf(end[counter]) + 7;
    rnaConverted =
      rnaConverted.slice(0, intronStart) + rnaConverted.slice(intronEnd);
    counter++;
  } while (counter < start.length);

  return rnaConverted;
};

const transcription = (dna) => {
  let rna = rnaPolymerase(dna);
  let mRna = splicing(rna);
  return mRna;
};

const translation = (mRna) => {
  let peptide = [];
  let translating = false;

  for (let i = 0; i < mRna.length; i += translating ? 3 : 1) {
    let codon = mRna.slice(i, i + 3);
    if (!translating) {
      if (codon === START_CODON) {
        translating = true;
        continue;
      } else {
        continue;
      }
    }
    let aminoAcid = CODONS[codon];
    if (aminoAcid === STOP_CODON) break;
    peptide.push(aminoAcid.key);
  }

  return peptide.join("");
};

const geneExpression = (gene) => {
  let count = 0;
  let mRna = transcription(gene);
  let protein = translation(mRna).split("");
  return protein.reduce((x, y) => {
    count++;
    if (count == 1 || count == protein.length || count == protein.length - 1) {
      return x + "";
    } else {
      return x + y;
    }
  }, "");
};

const geneExpressionMutated = (gene) => {
  let mutation = sickleCellAnemia(gene);
  gene = mutation;
  let count = 0;
  let mRna = transcription(gene);
  let protein = translation(mRna).split("");
  return protein.reduce((x, y) => {
    count++;
    if (count == 1 || count == protein.length || count == protein.length - 1) {
      return x + "";
    } else {
      return x + y;
    }
  }, "");
};

//------------------------------------------------------

const normalDisplay = document.querySelector(".expression__menu--normal");
const mutateDisplay = document.querySelector(".expression__menu--mutate");

normalDisplay.addEventListener("click", () => {
  const haemoglobinSubunitBeta = geneExpression(HBB_GENE);
  console.log(haemoglobinSubunitBeta);
});

mutateDisplay.addEventListener("click", () => {
  const haemoglobinSubunitBeta = geneExpressionMutated(HBB_GENE);
  console.log(haemoglobinSubunitBeta);
});
