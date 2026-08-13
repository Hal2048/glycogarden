const assert = require('node:assert/strict')

global.window = { location: { hostname: 'localhost' } }
global.Vue = {
  createApp: () => ({ mount: () => {} }),
  ref: value => ({ value }),
  computed: fn => ({ get value() { return fn() } }),
  onMounted: () => {},
  watch: () => {},
  nextTick: async () => {},
}

require('../app.js')

const { normalizeProfile, rebalanceProfile } = global.DistributionUtils
const compartments = ['CGC', 'MGC', 'TGC', 'TGN']

function almostEqual(actual, expected, tolerance = 1e-12) {
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`)
}

function assertConserved(profile) {
  almostEqual(compartments.reduce((sum, name) => sum + profile[name], 0), 1)
  compartments.forEach(name => assert.ok(profile[name] >= 0))
}

{
  const profile = { CGC: 0.1, MGC: 0.2, TGC: 0.3, TGN: 0.4 }
  const locks = { CGC: false, MGC: false, TGC: false, TGN: false }
  const result = rebalanceProfile(profile, 'CGC', 0.4, locks, compartments)
  almostEqual(result.CGC, 0.4)
  almostEqual(result.MGC, 0.6 * (0.2 / 0.9))
  almostEqual(result.TGC, 0.6 * (0.3 / 0.9))
  almostEqual(result.TGN, 0.6 * (0.4 / 0.9))
  assertConserved(result)
}

{
  const profile = { CGC: 0.1, MGC: 0.2, TGC: 0.3, TGN: 0.4 }
  const locks = { CGC: false, MGC: true, TGC: false, TGN: false }
  const result = rebalanceProfile(profile, 'CGC', 0.5, locks, compartments)
  almostEqual(result.MGC, 0.2)
  almostEqual(result.CGC, 0.5)
  assertConserved(result)
}

{
  const profile = { CGC: 1, MGC: 0, TGC: 0, TGN: 0 }
  const locks = { CGC: false, MGC: false, TGC: false, TGN: false }
  const result = rebalanceProfile(profile, 'CGC', 0.25, locks, compartments)
  almostEqual(result.MGC, 0.25)
  almostEqual(result.TGC, 0.25)
  almostEqual(result.TGN, 0.25)
  assertConserved(result)
}

{
  const profile = { CGC: 0.1, MGC: 0.2, TGC: 0.3, TGN: 0.4 }
  const locks = { CGC: false, MGC: true, TGC: true, TGN: true }
  const result = rebalanceProfile(profile, 'CGC', 0.8, locks, compartments)
  almostEqual(result.CGC, 0.1)
  assertConserved(result)
}

{
  const result = normalizeProfile({ CGC: 0, MGC: 0, TGC: 0, TGN: 0 }, compartments)
  compartments.forEach(name => almostEqual(result[name], 0.25))
}

console.log('distribution frontend tests passed')
