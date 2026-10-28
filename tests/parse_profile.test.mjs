// Unit tests for parseProfile() in src/js/worker.js (run with: node --test)
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

// worker.js is a Web Worker script, so pull the two parsing functions out of it.
const source = readFileSync(fileURLToPath(new URL('../src/js/worker.js', import.meta.url)), 'utf8');
function extract(name) {
    const start = source.indexOf(`function ${name}(`);
    let depth = 0;
    for (let i = source.indexOf('{', start); i < source.length; i++) {
        if (source[i] === '{') depth++;
        else if (source[i] === '}' && --depth === 0) return source.slice(start, i + 1);
    }
}
const { parseProfile } = new Function(
    `${extract('detectDelimiter')}\n${extract('parseProfile')}\nreturn { parseProfile };`)();

test('an ST / ST_id column is not a locus; # columns are skipped too', () => {
    const r = parseProfile('#Name\tST\tg1\tg2\n#ignored\ts1\t7\t1\t1\ns2\t9\t1\t2\n');
    assert.deepEqual(r.geneNames, ['g1', 'g2']);
    assert.equal(r.profiles[0].length, 2);
    const withStId = parseProfile('#Name\tst_Id\t#note\tg1\ns1\t5\tx\t1\n');
    assert.deepEqual(withStId.geneNames, ['g1']);
});

test('only 0, N and - (any case) are missing; -1 is an allele', () => {
    const r = parseProfile('#Name\tg1\tg2\tg3\tg4\tg5\ns1\t0\tn\t-\t-1\t1\ns2\t1\t1\t1\t1\t1\n');
    assert.deepEqual(r.profiles[0].slice(0, 3), [0, 0, 0]);
    assert.ok(r.profiles[0][3] > 0, '-1 must not be treated as missing');
    assert.notEqual(r.profiles[0][3], r.profiles[1][3], '-1 and 1 are different alleles');
    assert.equal(r.profiles[0][4], r.profiles[1][4]);
});

test('equal labels in a locus get equal ids, different labels different ids', () => {
    const r = parseProfile('#Name\tg1\ns1\t12\ns2\t12\ns3\t7\ns4\tabc\ns5\tABC\n');
    const p = r.profiles.map(row => row[0]);
    assert.equal(p[0], p[1]);
    assert.notEqual(p[0], p[2]);
    assert.equal(p[3], p[4]);       // labels are compared case-insensitively, as in GrapeTree
    assert.ok(p.every(v => v > 0));
});

test('CRLF and CR line endings, comma delimiter, header without #', () => {
    const crlf = parseProfile('#Name\tg1\tg2\r\ns1\t1\t2\r\ns2\t1\t3\r\n');
    assert.deepEqual(crlf.strains, ['s1', 's2']);
    const comma = parseProfile('Name,g1,g2\ns1,1,2\ns2,2,2\n');
    assert.deepEqual(comma.geneNames, ['g1', 'g2']);
    assert.equal(comma.profiles.length, 2);
});

test('empty and N/A fields are missing; malformed rows are skipped', () => {
    const r = parseProfile('#Name\tg1\tg2\ns1\t\tN/A\ns2\t1\ns3\t1\t1\n');
    assert.deepEqual(r.strains, ['s1', 's3']);
    assert.deepEqual(r.profiles[0], [0, 0]);
});

test('a header without any locus column is an error', () => {
    assert.throws(() => parseProfile('#Name\tST\ns1\t1\n'), /No locus columns/);
});
