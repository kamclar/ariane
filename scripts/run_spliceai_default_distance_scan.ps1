<# Run the offline D=50/500/4999 extension of a completed high-distance scan. #>
[CmdletBinding()]
param(
    [ValidateSet('Plan', 'Run', 'Resume', 'Status')]
    [string]$Action = 'Plan',
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$')]
    [string]$BaseOutputName = 'distance-stratified',
    [ValidatePattern('^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$')]
    [string]$OutputName = 'distance-defaults',
    [ValidateRange(0, 10000000)]
    [int]$MaxCases = 0,
    [ValidateRange(0, 604800)]
    [int]$MaxSeconds = 0,
    [switch]$Detached,
    [switch]$RetryErrors
)
$ErrorActionPreference = 'Stop'
$scanRoot = Split-Path -Parent $PSScriptRoot
$scanBuild = Join-Path $scanRoot 'data/spliceai/build'
$scanOutput = Join-Path $scanBuild $OutputName
if ($Action -eq 'Status') {
    $summaryPath = Join-Path $scanOutput 'summary.json'
    if (Test-Path -LiteralPath $summaryPath) {
        $summary = Get-Content -LiteralPath $summaryPath -Raw | ConvertFrom-Json
        $summary | Select-Object status, updated_at, counts, flags | ConvertTo-Json -Depth 5
    } else {
        Write-Output "No execution summary yet. Plan: $(Join-Path $scanOutput 'plan.json')"
    }
    exit 0
}
if ($RetryErrors -and $Action -ne 'Resume') { throw 'RetryErrors requires Resume.' }
if ($Detached -and $Action -eq 'Plan') { throw 'Detached is only for Run/Resume.' }
$profilePath = Join-Path $scanRoot 'data/spliceai/enigma_v1_2_spliceai_profile.json'
$scanImage = (Get-Content -LiteralPath $profilePath -Raw | ConvertFrom-Json).approved_engine.docker_image
New-Item -ItemType Directory -Path $scanBuild -Force | Out-Null
$dockerArgs = @('run', '--rm', '--init', '--network', 'none', '--cap-drop', 'ALL',
    '--security-opt', 'no-new-privileges:true', '--cpus', '4', '--memory', '6g',
    '--env', 'PYTHONDONTWRITEBYTECODE=1', '--env', 'OMP_NUM_THREADS=4',
    '--env', 'TF_NUM_INTRAOP_THREADS=4', '--env', 'TF_NUM_INTEROP_THREADS=1',
    '--mount', "type=bind,src=$scanRoot,dst=/workspace,readonly",
    '--mount', "type=bind,src=$scanBuild,dst=/experiment-output",
    '--workdir', '/workspace')
if ($Detached) { $dockerArgs += @('--detach', '--name', "ariane-$OutputName") }
$dockerArgs += @('--entrypoint', 'python', $scanImage,
    '-m', 'scripts.spliceai_default_distance_scan', '--fasta', '/hg38.fa.gz',
    '--base-output', "/experiment-output/$BaseOutputName",
    '--output', "/experiment-output/$OutputName", '--image-digest', $scanImage)
if ($Action -eq 'Plan') { $dockerArgs += '--prepare-only' }
if ($Action -eq 'Resume') { $dockerArgs += '--resume' }
if ($RetryErrors) { $dockerArgs += '--retry-errors' }
if ($MaxCases -gt 0) { $dockerArgs += @('--max-cases', "$MaxCases") }
if ($MaxSeconds -gt 0) { $dockerArgs += @('--max-seconds', "$MaxSeconds") }
& docker @dockerArgs
if ($LASTEXITCODE -ne 0) { throw "Docker experiment failed (exit $LASTEXITCODE)." }
