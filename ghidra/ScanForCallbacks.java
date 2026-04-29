import java.util.LinkedList;
import java.util.List;

import ghidra.app.script.GhidraScript;
import ghidra.program.model.address.Address;
import ghidra.program.model.address.AddressIterator;
import ghidra.program.model.listing.FunctionManager;
import ghidra.program.model.listing.Program;
import ghidra.program.model.symbol.ReferenceManager;

public class ScanForCallbacks extends GhidraScript{
	
	private Program program;
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		Address start_addr = this.askAddress("Start Address", "Set scan start address");
		Address end_addr = this.askAddress("End Address", "Set scan end address");
		
		ReferenceManager refMgr = program.getReferenceManager();
		FunctionManager funcMgr = program.getFunctionManager();
		
		AddressIterator allrefs = refMgr.getReferenceDestinationIterator(start_addr, true);
		List<Address> allAddr = new LinkedList<Address>(); //Doing this just in case adding functions messes up the iterator
		for(Address testAddr : allrefs) allAddr.add(testAddr);
		
		for(Address testAddr : allAddr) {
			if(testAddr.compareTo(end_addr) >= 0) break;
			
			if(funcMgr.getFunctionContaining(testAddr) == null) {
				//Reference or jump target that is not in a function
				//program.getMemory()
				if(getInstructionAt(testAddr) != null) {
					this.createFunction(testAddr, String.format("cbfunc_%08X", testAddr.getOffset()));
				}
			}
		}
		
		
		
	}

}
