/* BIO 116 genome browser lab: graded exercise engine, shared by index.html and grader.html.
   generate(code) builds a student's personal question set from the sequences in loci.js. The code is a
   random question-set code the student page makes up for each student; the same code always gives the
   same questions, so the grader can rebuild them from a submission file.
   score(submission) grades the auto-scored parts. The student page never calls score(). */
window.GBExercise = (() => {
'use strict';
const VERSION = 1;
const RAW = JSON.parse(JSON.stringify({
  hras: {chr: window.GB_LOCI.hras.chr, seq: window.GB_LOCI.hras.seq, genes: window.GB_LOCI.hras.genes},
  hbb:  {chr: window.GB_LOCI.hbb.chr,  seq: window.GB_LOCI.hbb.seq,  genes: window.GB_LOCI.hbb.genes}
}));
const COMP = {A:'T',T:'A',G:'C',C:'G'};
const BASES = 'TCAG', AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const tr = c => AAS[BASES.indexOf(c[0])*16 + BASES.indexOf(c[1])*4 + BASES.indexOf(c[2])];
const AA3 = {A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',Q:'Gln',E:'Glu',G:'Gly',H:'His',I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',T:'Thr',W:'Trp',Y:'Tyr',V:'Val','*':'Stop'};
const AAN = {A:'Alanine',R:'Arginine',N:'Asparagine',D:'Aspartate',C:'Cysteine',Q:'Glutamine',E:'Glutamate',G:'Glycine',H:'Histidine',I:'Isoleucine',L:'Leucine',K:'Lysine',M:'Methionine',F:'Phenylalanine',P:'Proline',S:'Serine',T:'Threonine',W:'Tryptophan',Y:'Tyrosine',V:'Valine','*':'Stop'};
const revcomp = s => s.split('').reverse().map(b => COMP[b]).join('');
const fmt = n => n.toLocaleString('en-US');

// Plus-strand coding genes only (every gene used here is + in its record).
function codingTx(locus, name){
  const L = RAW[locus], g = L.genes.find(g => g.name===name), t = g.txs[0];
  const cds = [];
  for (const [s,e] of t.exons) for (let p=s;p<=e;p++) if (p>=t.thick[0] && p<=t.thick[1]) cds.push(p);
  const base = p => L.seq[p-1];
  return {locus, chr:L.chr, name, id:t.id, exons:t.exons, cds, base, codon: k => cds.slice((k-1)*3, k*3).map(base).join('')};
}

// hg38 facts used in the UCSC question: chromosome and strand of each gene's MANE/RefSeq transcript.
const UCSC_GENES = [
  ['TP53','17','-'],['BRCA1','17','-'],['BRCA2','13','+'],['CFTR','7','+'],['KRAS','12','-'],
  ['EGFR','7','+'],['INS','11','-'],['DMD','X','-'],['APOE','19','+'],['MYC','8','+'],
  ['GAPDH','12','+'],['ACTB','7','-'],['HTT','4','+'],['SHH','7','-'],['ALB','4','+'],
  ['SOD1','21','+'],['APP','21','-'],['HBA1','16','+'],['LCT','2','-'],['FBN1','15','-']
];
const WRITTEN = [
  'IVS1-110 G>A sits 110 bases inside intron 1 of HBB, yet it causes β-thalassemia. Explain how, using what you saw in the browser.',
  'The browser calls the sickle-cell change codon 7 of HBB, while textbooks call it β6. Explain why, and say which numbering HGVS uses.',
  'In U01317 the β-globin genes run left to right, but in UCSC (hg38) HBB is the leftmost gene with arrows pointing left. Explain why.',
  'HRAS and LRRC56 share one CpG island in NG_007666. Describe how the two genes are arranged and why a CpG island often sits where it does.'
];

function normId(id){ return String(id||'').trim().toUpperCase().replace(/\s+/g,''); }
// Question-set codes: 6 characters, no look-alikes (0/O, 1/I/L).
const CODE_ABC = 'ABCDEFGHJKMNPQRSTUVWXYZ23456789';
function newCode(){
  const a = new Uint32Array(6);
  if (window.crypto?.getRandomValues) crypto.getRandomValues(a); else for (let i=0;i<6;i++) a[i] = Math.floor(Math.random()*2**32);
  return [...a].map(n => CODE_ABC[n % CODE_ABC.length]).join('');
}
const codeOf = sub => sub.setCode ?? sub.studentId;
function hash32(str){ let h = 0x811c9dc5; for (let i=0;i<str.length;i++){ h ^= str.charCodeAt(i); h = Math.imul(h, 0x01000193); } return h >>> 0; }
function rng(seed){ let a = seed|0; return () => { a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
const pick = (r, arr) => arr[Math.floor(r()*arr.length)];
const between = (r, lo, hi) => lo + Math.floor(r()*(hi-lo+1));

function generate(code){
  const id = normId(code);
  const r = rng(hash32('bio116-gb-v' + VERSION + ':' + id));
  const Q = [];

  // Q1 + Q2: a coding position in a globin gene, then a hypothetical change there
  const gname = pick(r, ['HBE1','HBG2','HBG1','HBD','HBB']);
  const gt = codingTx('hbb', gname), nCod = gt.cds.length/3;
  // codons 2-30 stay in the first coding exon (92 bp in every globin)
  const k = between(r, 2, Math.min(30, nCod-1)), j = between(r, 0, 2), p = gt.cds[(k-1)*3 + j];
  const ref = gt.base(p), refCodon = gt.codon(k);
  Q.push({id:'q1', title:'Read a coding position', locus:'hbb',
    prompt:`In the β-globin cluster tab, go to <code>${gt.chr}:${fmt(p)}</code> and zoom to base level.`,
    parts:[
      {id:'gene', label:'Which gene’s coding sequence is this base in?', kind:'text', key:gname},
      {id:'base', label:'Reference base at this position (+ strand)', kind:'base', key:ref},
      {id:'codon', label:'Codon number, counting the ATG start codon as 1', kind:'int', key:k},
      {id:'aa', label:'Amino acid encoded by that codon', kind:'aa', key:tr(refCodon)}
    ]});
  const alts = ['A','C','G','T'].filter(b => b!==ref), alt = pick(r, alts);
  const altCodon = refCodon.slice(0,j) + alt + refCodon.slice(j+1), aR = tr(refCodon), aA = tr(altCodon);
  const eff = aR===aA ? 'Synonymous' : aA==='*' ? 'Nonsense (stop gained)' : 'Missense';
  Q.push({id:'q2', title:'Predict a variant’s effect', locus:'hbb',
    prompt:`Suppose a patient has <b>${alt}</b> instead of <b>${ref}</b> at <code>${gt.chr}:${fmt(p)}</code> (the base from question 1).`,
    parts:[
      {id:'newcodon', label:'The new codon (3 letters, 5′→3′ on the gene’s coding strand)', kind:'dna3', key:altCodon},
      {id:'effect', label:'Effect on the protein', kind:'choice', opts:['Synonymous','Missense','Nonsense (stop gained)'], key:eff}
    ]});

  // Q3: an HRAS splice site
  const h = codingTx('hras', 'HRAS'), ex = h.exons, ik = between(r, 1, ex.length-1);
  const donor = r() < 0.5;
  const q3pos = donor ? ex[ik-1][1] + 1 : ex[ik][0] - 1;
  const q3bases = donor ? h.base(q3pos) + h.base(q3pos+1) : h.base(q3pos-1) + h.base(q3pos);
  Q.push({id:'q3', title:'Find a splice site', locus:'hras',
    prompt:`In the HRAS tab, find intron ${ik} of HRAS (NM_005343.4), between exon ${ik} and exon ${ik+1}.`,
    parts:[
      {id:'pos', label: donor ? `Coordinate of the FIRST base of intron ${ik}` : `Coordinate of the LAST base of intron ${ik}`, kind:'int', key:q3pos},
      {id:'bases', label: donor ? 'The first two bases of the intron' : 'The last two bases of the intron', kind:'dna', key:q3bases}
    ]});

  // Q4: genome orientation of an HRAS codon
  // exon 2 holds codons 1-37
  const hk = between(r, 2, 37), hc = h.codon(hk);
  Q.push({id:'q4', title:'Switch to genome orientation', locus:'hras',
    prompt:`NG_007666 shows HRAS on its coding strand, but in GRCh38 HRAS lies on the − strand of chr11. Find codon ${hk} of HRAS in the HRAS tab.`,
    parts:[
      {id:'codon', label:`Codon ${hk} as it reads in NG_007666 (5′→3′)`, kind:'dna3', key:hc},
      {id:'plus', label:`The same three bases read left to right on chr11’s + strand`, kind:'dna3', key:revcomp(hc)}
    ]});

  // Q5: BED coordinates of an exon
  const bLocus = r() < 0.5 ? 'hbb' : 'hras';
  const bName = bLocus==='hbb' ? pick(r, ['HBE1','HBG2','HBG1','HBD','HBB']) : 'HRAS';
  const bt = codingTx(bLocus, bName), en = between(r, 1, bt.exons.length), [es, ee] = bt.exons[en-1];
  Q.push({id:'q5', title:'Write a BED line', locus:bLocus,
    prompt:`Find exon ${en} of ${bName} in the ${bLocus==='hbb'?'β-globin cluster':'HRAS'} tab (click the gene; exons are listed 5′→3′).`,
    parts:[
      {id:'start', label:'BED start (0-based)', kind:'int', key:es-1},
      {id:'end', label:'BED end', kind:'int', key:ee},
      {id:'len', label:'Exon length in bp', kind:'int', key:ee-es+1}
    ]});

  // Q6: UCSC on the real genome
  const [ug, uchr, ustr] = pick(r, UCSC_GENES);
  Q.push({id:'q6', title:'Use the real genome (UCSC, hg38)', locus:null,
    prompt:`Open <a href="https://genome.ucsc.edu" target="_blank" rel="noopener">genome.ucsc.edu</a>, choose Human GRCh38/hg38, and search for <b>${ug}</b>.`,
    parts:[
      {id:'chr', label:'Chromosome', kind:'chr', key:uchr},
      {id:'strand', label:'Strand of the gene', kind:'choice', opts:['+','−'], key: ustr==='+' ? '+' : '−'},
      {id:'start', label:'At the start codon, which three letters does UCSC’s sequence track show, left to right?', kind:'dna3', key: ustr==='+' ? 'ATG' : 'CAT'},
      {id:'pos', label:'The position box after you zoom to the whole gene (copy it by typing; not auto-scored)', kind:'free', key:null}
    ]});

  // Q7: written, instructor-scored
  Q.push({id:'q7', title:'Explain it', locus:null, prompt: pick(r, WRITTEN),
    parts:[{id:'text', label:'Your answer (3–5 sentences)', kind:'long', key:null}]});

  return {version:VERSION, setCode:id, questions:Q};
}

/* ---------- grading ---------- */
const clean = s => String(s ?? '').trim();
const AA_NAMES = {};
for (const k of Object.keys(AA3)){ AA_NAMES[k.toUpperCase()] = k; AA_NAMES[AA3[k].toUpperCase()] = k; AA_NAMES[AAN[k].toUpperCase()] = k; }
Object.assign(AA_NAMES, {GLUTAMIC:'E', 'GLUTAMIC ACID':'E', 'ASPARTIC ACID':'D', STOP:'*', TER:'*', 'X':'*', 'STOP CODON':'*'});
function check(part, given){
  const g = clean(given);
  switch (part.kind){
    case 'text': return g.toUpperCase().replace(/\s+/g,'') === String(part.key).toUpperCase();
    case 'base': return g.toUpperCase() === part.key;
    case 'int':  return Number(g.replace(/[,\s]/g,'')) === part.key;
    case 'aa':   return AA_NAMES[g.toUpperCase().replace(/[()]/g,'').trim()] === part.key;
    case 'dna3': case 'dna': return g.toUpperCase().replace(/[^ACGTU]/g,'').replace(/U/g,'T') === part.key;
    case 'choice': return g === part.key || (part.key==='−' && g==='-');
    case 'chr':  return g.toLowerCase().replace(/^chr/,'').toUpperCase() === part.key;
    default: return null; // not auto-scored
  }
}
function keyText(part){
  if (part.key==null) return '(instructor)';
  if (part.kind==='aa') return `${AAN[part.key]} (${AA3[part.key]}, ${part.key})`;
  if (part.kind==='int') return fmt(part.key);
  return String(part.key);
}
// A light integrity check: flags files edited by hand after they were saved.
function seal(sub){
  const body = JSON.stringify([sub.version, codeOf(sub), sub.name, sub.answers, sub.startedAt, sub.submittedAt, sub.pasteAttempts]);
  return hash32('bio116-seal:' + body).toString(36) + hash32(body + ':bio116').toString(36);
}
function score(sub){
  const set = generate(codeOf(sub)), rows = [];
  let got = 0, max = 0;
  for (const q of set.questions) for (const part of q.parts){
    const given = (sub.answers?.[q.id] || {})[part.id];
    const ok = check(part, given);
    if (ok !== null){ max++; if (ok) got++; }
    rows.push({q:q.id, title:q.title, part:part.id, label:part.label, given: clean(given), key:keyText(part), ok});
  }
  return {got, max, rows, sealOk: sub.seal === seal(sub), versionOk: sub.version === VERSION};
}
return {VERSION, generate, score, seal, normId, newCode, codeOf};
})();
