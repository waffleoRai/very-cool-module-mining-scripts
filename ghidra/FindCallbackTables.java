import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import ghidra.app.script.GhidraScript;
import ghidra.program.model.address.Address;
import ghidra.program.model.address.AddressIterator;
import ghidra.program.model.data.DataType;
import ghidra.program.model.listing.Data;
import ghidra.program.model.listing.Program;
import ghidra.program.model.mem.Memory;
import ghidra.program.model.symbol.ReferenceManager;

public class FindCallbackTables extends GhidraScript{
	
	private Program program;
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		File outfile = this.askFile("Select Output File (csv)", "Okay");
		Address stAddr = this.askAddress("Start Address", "Enter start of data section to scan");
		Address edAddr = this.askAddress("End Address", "Enter end of data section to scan");
		
		ReferenceManager refMgr = program.getReferenceManager();
		Memory memMgr = program.getMemory();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("ADDRESS,ELEMENT_COUNT,CONTENTS\n"); //Contents are addresses, not symbols
		//Members are divided by ;
		AddressIterator allrefs = refMgr.getReferenceDestinationIterator(stAddr, true);
		for(Address testAddr : allrefs) {
			if(testAddr.compareTo(edAddr) >= 0) break;

			int symSize = 0;
			Address otherAddr = testAddr.add(1);
			while((otherAddr.compareTo(edAddr) < 0) && (!refMgr.hasReferencesTo(otherAddr))) {
				symSize++;
				otherAddr = otherAddr.add(1);
			}
			
			if(symSize < 4) continue;
			
			long addrOff = testAddr.getOffset();
			Data dd = this.getDataAt(testAddr);
			if(dd == null) continue;
			DataType dt = dd.getDataType();
			if(dt == null) continue;
			if(dt.getLength() != 4) continue;
			if((addrOff & 0x3L) == 0) {
				int maxWords = symSize >>> 2;
				int pointerCount = 0;
				int funcPtrCount = 0;
				int[] dat = new int[maxWords];
				memMgr.getInts(testAddr, dat);
				for(int i = 0; i < dat.length; i++) {
					if(dat[i] != 0) {
						Address trgAddr = this.toAddr(dat[i]);
						if(trgAddr == null) break;
						if(this.getFunctionAt(trgAddr) == null) break;
						funcPtrCount++;
					}
					pointerCount = i + 1;
				}
				
				if((funcPtrCount > 0) && (pointerCount > 1)) {
					bw.write(String.format("0x%08x,", testAddr.getOffset()));
					bw.write(pointerCount + ",");
					for(int i = 0; i < dat.length; i++) {
						if(i > 0) bw.write(";");
						bw.write(String.format("0x%08x,", dat[i]));
					}
					bw.write("\n");
				}
			}
		}
		
		bw.close();
	}

}
