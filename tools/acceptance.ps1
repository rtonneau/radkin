<#
.SYNOPSIS
  Phase-1 acceptance matrix: radkin vs radiopyo on the bundled configs, every
  beam type, single and chunked runs (tools/compare_radiopyo.py for each case).

.EXAMPLE
  pwsh tools/acceptance.ps1 -Python E:/miniconda3/envs/radiobio/python.exe -Radiopyo D:/DEV/PYTHON/radiopyo
#>
param(
    [string]$Python = $(if ($env:RADKIN_PYTHON) { $env:RADKIN_PYTHON } else { "python" }),
    [string]$Radiopyo = $(if ($env:RADIOPYO_DIR) { $env:RADIOPYO_DIR } else { "../../PYTHON/radiopyo" }),
    [double]$Rtol = 1e-6
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Push-Location $root
try {
    cargo build --release --quiet
    if ($LASTEXITCODE -ne 0) { throw "cargo build failed" }

    $configs = @(
        @{ Name = "pure_water"; Path = "data/advanced_water_radiolysis.toml"; Extra = @() },
        @{ Name = "labarbe"; Path = "data/config_Labarbe.toml"; Extra = @() },
        @{ Name = "paper"; Path = "data/CONFIG_PAPER_01.toml"; Extra = @("--variable-gvalues", "--let", "0.56") }
    )
    # Each beam: parameters, t_span, and per mode the step options.
    $beams = @(
        @{ Type = "constant"; Params = @("dose_rate=1e3", "max_dose=10"); Span = @("0", "1");
           Run = @(); Chunked = @("--max-step-on", "1e-3", "--max-step-off", "0.05") },
        @{ Type = "singlepulse"; Params = @("dose_rate=1e6", "max_dose=10"); Span = @("0", "1");
           Run = @("--max-step", "1e-6"); Chunked = @("--max-step-on", "1e-7", "--max-step-off", "0.05") },
        @{ Type = "pulsed"; Params = @("dose_rate_peak=1e6", "period=0.01", "on_time=2e-6", "max_dose=10"); Span = @("0", "0.1");
           RunSpan = @("0", "0.015"); Run = @("--max-step", "1e-6"); Chunked = @("--max-step-on", "1e-7", "--max-step-off", "1e-3") },
        @{ Type = "macro_pulsed"; Params = @("dose_rate_peak=1e6", "period=0.05", "on_time=0.001", "micro_period=1e-4", "micro_on_time=2e-6", "max_dose=10"); Span = @("0", "0.1");
           RunSpan = @("0", "0.0015"); Run = @("--max-step", "1e-6"); Chunked = @("--max-step-on", "1e-7", "--max-step-off", "1e-4") }
    )
    $results = @()
    foreach ($c in $configs) {
        foreach ($b in $beams) {
            $modes = @("run", "chunked")
            if ($b.Type -eq "singlepulse") { $modes += "single-pulse" }
            foreach ($mode in $modes) {
                $beamArgs = @("--beam-type", $b.Type)
                foreach ($p in $b.Params) { $beamArgs += @("--beam-param", $p) }
                if ($mode -eq "run") {
                    $span = if ($b.RunSpan) { $b.RunSpan } else { $b.Span }
                    $opts = $b.Run
                } elseif ($mode -eq "single-pulse") {
                    $span = $b.Span
                    $opts = @("--max-step-on", "1e-7", "--max-step-trans", "1e-5", "--max-step-final", "0.05")
                } else {
                    $span = $b.Span
                    $opts = $b.Chunked
                }
                $cmd = @("-I", "tools/compare_radiopyo.py", $c.Path, "--mode", $mode, "--t-span") + $span +
                    $beamArgs + $opts + $c.Extra + @("--rtol", "$Rtol", "--radiopyo", $Radiopyo)
                $out = & $Python @cmd 2>&1 | Where-Object { $_ -notmatch "INFO" }
                $ok = $LASTEXITCODE -eq 0
                $label = "{0,-10} {1,-12} {2,-12}" -f $c.Name, $b.Type, $mode
                Write-Host ("{0} {1}" -f $label, $(if ($ok) { "PASS" } else { "FAIL" }))
                if (-not $ok) { $out | ForEach-Object { Write-Host "    $_" } }
                $results += [pscustomobject]@{ Case = $label; Pass = $ok }
            }
        }
    }
    $failed = @($results | Where-Object { -not $_.Pass })
    Write-Host ""
    Write-Host ("{0}/{1} cases passed" -f ($results.Count - $failed.Count), $results.Count)
    if ($failed.Count -gt 0) { exit 1 }
}
finally {
    Pop-Location
}
