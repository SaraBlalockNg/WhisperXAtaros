# WhisperXAtaros
Utility scripts for processing the output of Whisper ASR and comparing to ATAROS-style textgrids. These scripts assume that:
1) ATAROS-style TextGrids are stored together in a single folder, where the first speaker ID NW{F|M}{000-999} is recorded first in the TextGrid
2) Outputs from running Whisper ASR on the audio are stored together in a single folder. In particular, that the `score.txt` files are stored in this folder.
3) For processing audio, that `sox` is an available command.

First, download required packages:
`pip install -r requirements.txt`

To align all transcripts in a set of folders, run:
`python align_ataros_candidates.py --indir {PATH TO WHISPER} --tg_path {PATH TO ATAROS}`
