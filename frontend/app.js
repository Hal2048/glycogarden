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

function cloneDistribution(distribution) {
  return JSON.parse(JSON.stringify(distribution || {}))
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
    const selectedPromoter = ref('')
    const selectedPreset = ref('baseline')
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
            return `${row.label}<br/>Rank ${row.rank}: ${formatPercent(p.value)}`
          },
        },
        xAxis: {
          type: 'category',
          data: rows.map(d => d.label),
          axisLabel: { rotate: 45, interval: 0, fontSize: 10 },
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
        selectedPromoter.value = config.value.promoterLevels[0]?.key || ''
        enzymeDistribution.value = cloneDistribution(config.value.baselineDistribution)
        distributionLocks.value = emptyLocks()
      } catch (e) {
        error.value = e.message
      } finally {
        loadingConfig.value = false
      }
    }

    function selectPromoter(key) {
      if (selectedPromoter.value !== key) {
        selectedPromoter.value = key
        markResultsStale()
      }
    }

    function applyPreset(presetKey) {
      const preset = config.value?.enzymePresets.find(item => item.key === presetKey)
      if (!preset) return
      enzymeDistribution.value = cloneDistribution(preset.distribution)
      distributionLocks.value = emptyLocks()
      selectedPreset.value = presetKey
      inputError.value = ''
      markResultsStale()
    }

    function resetAll() {
      applyPreset('baseline')
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
      selectedPreset.value = 'custom'
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
      selectedPreset.value = 'custom'
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

    async function runSimulation() {
      if (!selectedPromoter.value || computing.value) return
      computing.value = true
      error.value = ''
      inputError.value = ''
      cancelledMessage.value = ''
      activeController = new AbortController()
      const requestDistribution = cloneDistribution(enzymeDistribution.value)
      const requestSignature = JSON.stringify(requestDistribution)
      try {
        const res = await fetch(`${API_URL}/predict`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            promoter: selectedPromoter.value,
            enzymeDistribution: requestDistribution,
          }),
          signal: activeController.signal,
        })
        const data = await res.json()
        if (!res.ok) throw new Error(data.error || `Prediction failed: ${res.status}`)
        result.value = data
        const controlsChanged = requestSignature !== JSON.stringify(enzymeDistribution.value)
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

    return {
      config,
      selectedPromoter,
      selectedPreset,
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
      selectPromoter,
      applyPreset,
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
    }
  },
}).mount('#app')
