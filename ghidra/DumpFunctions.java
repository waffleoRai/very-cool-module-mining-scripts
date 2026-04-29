import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import ghidra.app.script.GhidraScript;
import ghidra.program.model.address.Address;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.FunctionIterator;
import ghidra.program.model.listing.Program;
import ghidra.program.model.symbol.SourceType;

public class DumpFunctions extends GhidraScript{
	
	private Program program;
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		File outfile = this.askFile("Select Output File", "Okay");
		Address start_addr = this.askAddress("Start Address", "Set scan start address");
		Address end_addr = this.askAddress("End Address", "Set scan end address");
		
		FunctionIterator funcs = getCurrentProgram().getFunctionManager().getFunctions(start_addr, true);
		if(funcs == null) {
			printerr("No functions found in range...");
			return;
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
		out.write("#NAME\tADDRESS\tCALLTYPE\tSIG\tSIZE\n");
		Function last_func = null;
		for(Function func : funcs) {
			Address fstart = func.getEntryPoint();
			if(fstart.compareTo(end_addr) > 0) break;
			
			if(last_func != null) {
				out.write("0x");
				out.write(Long.toHexString(fstart.subtract(last_func.getEntryPoint())));
				out.write('\n');
			}
			//Address fend = func.getBody().getMaxAddress();
			String fname = func.getName();
			if(fname.startsWith("FUN_")) {
				fname = "func_" + fstart.toString().toUpperCase();
				func.setName(fname, SourceType.DEFAULT);
			}
			out.write(fname); out.write('\t');
			out.write(fstart.toString().toUpperCase()); out.write('\t');
			out.write(func.getCallingConventionName()); out.write('\t');
			out.write(func.getSignature().getPrototypeString()); out.write('\t');
			
			last_func = func;
		}
		if(last_func != null) {
			out.write("0x");
			out.write(Long.toHexString(end_addr.subtract(last_func.getEntryPoint())));
			out.write('\n');
		}
		out.close();
	}

}
