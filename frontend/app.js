const { createApp, ref, computed, onMounted, watch, nextTick } = Vue

const API_URL = window.location.hostname === 'localhost' ? 'http://localhost:5000/api' : '/api'
const EPSILON = 1e-12
const COMPARTMENT_COLORS = {
  CGC: '#6f8e6f',
  MGC: '#9caf88',
  TGC: '#d4a96a',
  TGN: '#a97979',
}

function formatPercent(value) {
  return `${(value * 100).toFixed(2)}%`
}

function formatConcentration(value) {
  if (!Number.isFinite(value)) return '—'
  if (value >= 1000) return value.toFixed(1)
  if (value >= 1) return value.toFixed(2)
  return value.toPrecision(3)
}

function cloneDistribution(distribution) {
  return JSON.parse(JSON.stringify(distribution || {}))
}

function defaultParams(config) {
  return {
    tau: config.tauDefault,
    compartmentVolume: config.compartmentVolumeDefault,
    proteinProdRate: config.proteinProdRateDefault,
    donorConcs: { ...config.donorDefaults },
    enzymeConcs: { ...config.enzymeConcDefaults },
  }
}

function normalizeProfile(profile, compartments) {
  const values = compartments.map(name => Math.max(0, Number(profile[name]) || 0))
  const total = values.reduce((sum, value) => sum + value, 0)
  if (total <= EPSILON) {
    const share = 1 / compartments.length
    return Object.fromEntries(compartments.map(name => [name, share]))
  }
  return Object.fromEntries(compartments.map((name, index) => [name, values[index] / total]))
}

function rebalanceProfile(profile, changedCompartment, requestedValue, locks, compartments) {
  const next = { ...profile }
  const peers = compartments.filter(name => name !== changedCompartment)
  const lockedPeers = peers.filter(name => locks[name])
  const adjustablePeers = peers.filter(name => !locks[name])
  const lockedTotal = lockedPeers.reduce((sum, name) => sum + next[name], 0)

  if (adjustablePeers.length === 0) {
    next[changedCompartment] = Math.max(0, 1 - lockedTotal)
    return normalizeProfile(next, compartments)
  }

  const changedValue = Math.min(Math.max(requestedValue, 0), Math.max(0, 1 - lockedTotal))
  const remainder = Math.max(0, 1 - lockedTotal - changedValue)
  const adjustableTotal = adjustablePeers.reduce((sum, name) => sum + next[name], 0)

  next[changedCompartment] = changedValue
  if (adjustableTotal > EPSILON) {
    adjustablePeers.forEach(name => {
      next[name] = remainder * (next[name] / adjustableTotal)
    })
  } else {
    const share = remainder / adjustablePeers.length
    adjustablePeers.forEach(name => {
      next[name] = share
    })
  }

  const residual = 1 - compartments.reduce((sum, name) => sum + next[name], 0)
  next[adjustablePeers[adjustablePeers.length - 1]] += residual
  return next
}

globalThis.DistributionUtils = { normalizeProfile, rebalanceProfile }

function makeChartRows(result) {
  if (!result || !result.top) return []
  return result.top.map(item => {
    const structure = result.structures?.find(s => s.id === item.id)
    return {
      id: item.id,
      label: structure?.label || `Glycoform ${item.id}`,
      value: item.value,
      rank: item.rank,
    }
  })
}

createApp({
  setup() {
    const config = ref(null)
    const params = ref(null)
    const enzymeDistribution = ref({})
    const distributionLocks = ref({})
    const result = ref(null)
    const resultStale = ref(false)
    const loadingConfig = ref(true)
    const computing = ref(false)
    const error = ref('')
    const inputError = ref('')
    const cancelledMessage = ref('')
    const chartRef = ref(null)
    let chart = null
    let activeController = null

    const chartRows = computed(() => makeChartRows(result.value))

    const totGlycanConc = computed(() => {
      const p = params.value
      if (!p) return NaN
      return (p.proteinProdRate * p.tau) / p.compartmentVolume
    })

    function emptyLocks() {
      return Object.fromEntries(config.value.enzymeNames.map(enzyme => [
        enzyme,
        Object.fromEntries(config.value.compartments.map(compartment => [compartment, false])),
      ]))
    }

    function markResultsStale() {
      if (result.value) resultStale.value = true
    }

    function renderChart() {
      if (!chartRef.value) return
      if (!chart) chart = echarts.init(chartRef.value)
      const rows = chartRows.value
      if (rows.length === 0) return
      chart.setOption({
        grid: { top: 40, right: 24, bottom: 80, left: 64 },
        tooltip: {
          trigger: 'axis',
          axisPointer: { type: 'shadow' },
          formatter: params => {
            const p = Array.isArray(params) ? params[0] : params
            const row = rows[p.dataIndex]
            return `ID ${row.id}<br/>${row.label}<br/>Rank ${row.rank}: ${formatPercent(p.value)}`
          },
        },
        xAxis: {
          type: 'category',
          data: rows.map(d => `#${d.id}`),
          axisLabel: { interval: 0, fontSize: 11 },
        },
        yAxis: {
          type: 'value',
          name: 'Relative abundance',
          axisLabel: { formatter: v => `${(v * 100).toFixed(0)}%` },
        },
        series: [{
          type: 'bar',
          data: rows.map(d => ({
            value: d.value,
            itemStyle: { color: `rgba(111, 142, 111, ${0.5 + 0.5 * (1 - d.rank / (rows.length + 1))})` },
          })),
          barMaxWidth: 32,
          label: {
            show: true,
            position: 'top',
            formatter: p => `${(p.value * 100).toFixed(1)}%`,
            fontSize: 10,
          },
        }],
      }, true)
    }

    async function loadConfig() {
      try {
        const res = await fetch(`${API_URL}/config`)
        if (!res.ok) throw new Error(`Failed to load config: ${res.status}`)
        config.value = await res.json()
        params.value = defaultParams(config.value)
        enzymeDistribution.value = cloneDistribution(config.value.baselineDistribution)
        distributionLocks.value = emptyLocks()
      } catch (e) {
        error.value = e.message
      } finally {
        loadingConfig.value = false
      }
    }

    function resetAll() {
      enzymeDistribution.value = cloneDistribution(config.value.baselineDistribution)
      distributionLocks.value = emptyLocks()
      params.value = defaultParams(config.value)
      inputError.value = ''
      markResultsStale()
    }

    function resetEnzyme(enzymeName) {
      enzymeDistribution.value = {
        ...enzymeDistribution.value,
        [enzymeName]: cloneDistribution(config.value.baselineDistribution[enzymeName]),
      }
      distributionLocks.value = {
        ...distributionLocks.value,
        [enzymeName]: Object.fromEntries(config.value.compartments.map(name => [name, false])),
      }
      params.value = {
        ...params.value,
        enzymeConcs: {
          ...params.value.enzymeConcs,
          [enzymeName]: config.value.enzymeConcDefaults[enzymeName],
        },
      }
      inputError.value = ''
      markResultsStale()
    }

    function toggleLock(enzymeName, compartment) {
      distributionLocks.value = {
        ...distributionLocks.value,
        [enzymeName]: {
          ...distributionLocks.value[enzymeName],
          [compartment]: !distributionLocks.value[enzymeName][compartment],
        },
      }
    }

    function commitCompartment(enzymeName, compartment, event) {
      const rawValue = event.target.value
      const percentage = Number(rawValue)
      if (rawValue.trim() === '' || !Number.isFinite(percentage) || percentage < 0 || percentage > 100) {
        inputError.value = `${enzymeName} ${compartment} must be a number from 0% to 100%.`
        event.target.value = displayPercentage(enzymeName, compartment)
        return
      }

      const profile = enzymeDistribution.value[enzymeName]
      const nextProfile = rebalanceProfile(
        profile,
        compartment,
        percentage / 100,
        distributionLocks.value[enzymeName],
        config.value.compartments,
      )
      enzymeDistribution.value = { ...enzymeDistribution.value, [enzymeName]: nextProfile }
      inputError.value = ''
      markResultsStale()
    }

    function displayPercentage(enzymeName, compartment) {
      return (enzymeDistribution.value[enzymeName]?.[compartment] * 100).toFixed(2)
    }

    function profileTotal(enzymeName) {
      const profile = enzymeDistribution.value[enzymeName]
      if (!profile) return '0.00%'
      const total = config.value.compartments.reduce((sum, name) => sum + profile[name], 0)
      return `${(total * 100).toFixed(2)}%`
    }

    function segmentStyle(enzymeName, compartment) {
      return {
        width: `${enzymeDistribution.value[enzymeName][compartment] * 100}%`,
        backgroundColor: COMPARTMENT_COLORS[compartment],
      }
    }

    function currentPayloadSignature() {
      return JSON.stringify({
        distribution: enzymeDistribution.value,
        params: params.value,
      })
    }

    async function runSimulation() {
      if (computing.value) return
      computing.value = true
      error.value = ''
      inputError.value = ''
      cancelledMessage.value = ''
      activeController = new AbortController()
      const requestDistribution = cloneDistribution(enzymeDistribution.value)
      const requestParams = JSON.parse(JSON.stringify(params.value))
      const requestSignature = currentPayloadSignature()
      try {
        const res = await fetch(`${API_URL}/predict`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            enzymeDistribution: requestDistribution,
            enzymeConcs: requestParams.enzymeConcs,
            tau: requestParams.tau,
            compartmentVolume: requestParams.compartmentVolume,
            proteinProdRate: requestParams.proteinProdRate,
            donorConcs: requestParams.donorConcs,
          }),
          signal: activeController.signal,
        })
        const data = await res.json()
        if (!res.ok) throw new Error(data.error || `Prediction failed: ${res.status}`)
        result.value = data
        const controlsChanged = requestSignature !== currentPayloadSignature()
        if (!controlsChanged) {
          enzymeDistribution.value = cloneDistribution(data.enzymeDistribution)
        }
        resultStale.value = controlsChanged
      } catch (e) {
        if (e.name === 'AbortError') {
          cancelledMessage.value = 'Stopped waiting for this run. The server may finish it and cache the result.'
        } else {
          error.value = e.message
        }
      } finally {
        computing.value = false
        activeController = null
        await nextTick()
        renderChart()
      }
    }

    function cancelSimulation() {
      activeController?.abort()
    }

    onMounted(() => {
      loadConfig()
      window.addEventListener('resize', () => chart?.resize())
    })

    watch([result, computing], async ([newResult, newComputing]) => {
      if (newResult && !newComputing) {
        await nextTick()
        renderChart()
      }
    })

    // Any edit to the physiology parameters invalidates the displayed result.
    watch(
      () => params.value && {
        tau: params.value.tau,
        v: params.value.compartmentVolume,
        q: params.value.proteinProdRate,
        d: params.value.donorConcs,
        e: params.value.enzymeConcs,
      },
      () => markResultsStale(),
    )

    return {
      config,
      params,
      enzymeDistribution,
      distributionLocks,
      result,
      resultStale,
      loadingConfig,
      computing,
      error,
      inputError,
      cancelledMessage,
      chartRef,
      chartRows,
      totGlycanConc,
      resetAll,
      resetEnzyme,
      toggleLock,
      commitCompartment,
      displayPercentage,
      profileTotal,
      segmentStyle,
      runSimulation,
      cancelSimulation,
      formatPercent,
      formatConcentration,
    }
  },
}).mount('#app')
