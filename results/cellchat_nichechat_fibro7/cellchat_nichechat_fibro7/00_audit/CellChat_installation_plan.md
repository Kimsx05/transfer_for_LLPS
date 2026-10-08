# Isolated CellChat installation record

State: COMPLETED AND LOAD-VERIFIED. Installation was authorized by the 2026-10-08 environment-installation request.

- Official source: https://github.com/jinworks/CellChat/releases/tag/v2.1.2
- Verified release commit: a0d3b2d231d46c8787177fffeac908270c253747.
- Installed version: 2.1.2.
- Library: /home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/Rlib.
- Missing dependencies ggalluvial, svglite, sna, ggnetwork were installed here. Existing project libraries were not upgraded.
- Installation log: CellChat_installation.log; fresh load verification and session: sessionInfo_CellChat_verified.txt.
- The installer log ends with an extraneous parse error after its PASS message because its comment was edited while the interpreter was reading the file. A fresh separate R process loaded the completed package and verified version 2.1.2 successfully. No incomplete installation was used for inference.
- NicheNet is separately installed in Rlib_nichenet. Its ggplot2 compatibility update does not enter the CellChat library path.
- NicheChat was a corrected name and is NOT_APPLICABLE; it is not a prerequisite or analysis method.
