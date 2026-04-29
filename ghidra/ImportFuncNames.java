import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;


import ghidra.app.script.GhidraScript;
import ghidra.program.model.address.Address;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.Program;
import ghidra.program.model.symbol.SourceType;

public class ImportFuncNames extends GhidraScript {

private Program program;
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		File infile = this.askFile("Select Input tsv", "Okay");
		//Expects NAME\tADDRESS
		int count = 0;
		BufferedReader br = new BufferedReader(new FileReader(infile));
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.isBlank()) continue;
			if(line.startsWith("#")) continue;
			String[] split = line.split("\t");
			if(split.length < 2) continue;
			
			String fname = split[0];
			String araw = split[1];
			long aval = 0L;
			
			try {
				if(araw.startsWith("0x")) araw = araw.substring(2);
				aval = Long.parseUnsignedLong(araw, 16);
			}
			catch(NumberFormatException ex) {
				ex.printStackTrace();
				continue;
			}
			
			Address myaddr = program.getAddressFactory().getAddress(araw);
			Function func = program.getFunctionManager().getFunctionAt(myaddr);
			if(func == null) {
				println("No function found at: " + myaddr.toString());
				continue;
			}
			func.setName(fname, SourceType.USER_DEFINED);
			
			count++;
		}
		br.close();
		
		println(count + " records processed!");
	}
	
}
