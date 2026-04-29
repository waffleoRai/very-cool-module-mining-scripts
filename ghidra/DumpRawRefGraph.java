import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import ghidra.app.script.GhidraScript;
import ghidra.program.model.address.Address;
import ghidra.program.model.address.AddressIterator;
import ghidra.program.model.data.DataType;
import ghidra.program.model.listing.Data;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.FunctionIterator;
import ghidra.program.model.listing.Program;
import ghidra.program.model.mem.Memory;
import ghidra.program.model.symbol.Reference;
import ghidra.program.model.symbol.ReferenceIterator;
import ghidra.program.model.symbol.ReferenceManager;
import ghidra.program.model.symbol.Symbol;
import ghidra.program.model.symbol.SymbolTable;

public class DumpRawRefGraph extends GhidraScript{
	
	private static class SectionDef{
		public String name;
		public Address start; //inclusive
		public Address end; //inclusive
		
		public boolean isInSection(Address addr) {
			if(addr.compareTo(start) < 0) return false;
			if(addr.compareTo(end) > 0) return false;
			return true;
		}
		
		public boolean isInSection(long addr) {
			if(addr < start.getOffset()) return false;
			if(addr > end.getOffset()) return false;
			return true;
		}
	}
	
	private Program program;
	
	private Map<String, SectionDef> readSecDefTable(File infile) throws IOException{
		Map<String, SectionDef> map = new HashMap<String, SectionDef>();
		
		//Expects tsv with three columns: Name Start End
		BufferedReader br = new BufferedReader(new FileReader(infile));
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.isEmpty()) continue;
			if(line.startsWith("#")) continue;
			
			String[] fields = line.split("\t");
			if(fields.length < 3) continue;
			SectionDef sec = new SectionDef();
			sec.name = fields[0].trim();
			String field = fields[1].trim();
			if(!field.startsWith("0x")) field = "0x" + field;
			sec.start = program.getAddressFactory().getAddress(field);
			field = fields[2].trim();
			if(!field.startsWith("0x")) field = "0x" + field;
			sec.end = program.getAddressFactory().getAddress(field);
			map.put(sec.name, sec);
		}
		br.close();
		return map;
	}
	
	private SectionDef getSection(Map<String, SectionDef> map, Address addr) {
		for(SectionDef def : map.values()) {
			if(def.isInSection(addr)) return def;
		}
		return null;
	}
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		File infile = this.askFile("Select Segment Def Table (tsv)", "Okay"); //Describes start/end of .text, .data, .bss etc.
		File outfile = this.askFile("Select Output File (tsv)", "Okay");
		
		Map<String, SectionDef> secTable = readSecDefTable(infile);
		SectionDef textDef = secTable.get(".text");
		
		if(textDef == null) {
			printerr(".text section was not defined in input table!");
			return;
		}
		
		//
		println("Scanning for unlabeled callback functions...");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("#SYMBOL\tADDRESS\tSECTION\tREFEREES\tTYPE\n"); //Referees are just addresses
		//bw.write("FromAddrHex,ToAddrHex,FromAddrSec,ToAddrSec\n");
		SymbolTable symTable = program.getSymbolTable();
		ReferenceManager refMgr = program.getReferenceManager();
		//Memory memMgr = program.getMemory();
		
		//.text
		FunctionIterator funcs = program.getFunctionManager().getFunctions(textDef.start, true);
		this.println("Working on .text...");
		for(Function func : funcs) {
			Address fstart = func.getEntryPoint();
			if(fstart.compareTo(textDef.end) >= 0) break;
			
			bw.write(func.getName());
			bw.write(String.format("\t0x%08x", fstart.getOffset()));
			bw.write("\t.text");
			
			bw.write("\t");
			Set<Long> checked = new HashSet<Long>();
			Reference[] refs = this.getReferencesTo(fstart);
			if((refs != null) && (refs.length > 0)) {
				int written = 0;
				checked.clear();
				for(int i = 0; i < refs.length; i++) {
					if(refs[i] == null) continue;
					Reference ref = refs[i];
					Address refereeAddr = ref.getFromAddress();
					long refOff = refereeAddr.getOffset();
					SectionDef refereeSec = getSection(secTable, refereeAddr);
					if(refereeSec != null && !checked.contains(refOff)) {
						if(written > 0) bw.write(";");
						bw.write(String.format("0x%08x", refereeAddr.getOffset()));
						written++;
					}
				}
			}
			bw.write("\t" + func.getCallingConventionName() + "\n");
		}
		
		//Other sections
		for(SectionDef def : secTable.values()) {
			if(def.name.equals(".text")) continue;
			
			this.println("Working on " + def.name + "...");
			
			//Try references directly?
			String anonNameStem = def.name;
			if(def.name.startsWith(".")) anonNameStem = def.name.substring(1);
			AddressIterator allrefs = refMgr.getReferenceDestinationIterator(def.start, true);
			for(Address testAddr : allrefs) {
				if(testAddr.compareTo(def.end) >= 0) break;
				
				Symbol[] syms = symTable.getSymbols(testAddr);
				Symbol sym = null;
				if((syms != null) && (syms.length > 0) && (syms[0] != null)) {
					sym = syms[0];
					bw.write(sym.getName());
				}
				else {
					bw.write(String.format("_%s_%08X", anonNameStem, testAddr.getOffset()));
				}
				
				bw.write(String.format("\t0x%08x", testAddr.getOffset()));
				bw.write("\t" + def.name);
				bw.write("\t");
				
				ReferenceIterator refs = refMgr.getReferencesTo(testAddr);
				if((refs != null) && refs.hasNext()) {
					int written = 0;
					for(Reference ref : refs) {
						Address refereeAddr = ref.getFromAddress();
						SectionDef refereeSec = getSection(secTable, refereeAddr);
						if(refereeSec != null) {
							if(written > 0) bw.write(";");
							bw.write(String.format("0x%08x", refereeAddr.getOffset()));
							written++;
						}
					}
				}
				
				//Try to fetch data type and array length
				bw.write("\t");
				Data dd = this.getDataAt(testAddr);
				if(dd != null) {
					DataType tt = dd.getDataType();
					if(tt != null) {
						bw.write("\t" + tt.getDisplayName());
						//TODO Maybe also check array len.
					}
					else bw.write("<UNK>");
				}
				else bw.write("<UNK>");
				bw.write("\n");
				
			}
			
		}

		bw.close();
	}
	

}
