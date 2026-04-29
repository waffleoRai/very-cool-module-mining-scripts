import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedList;
import java.util.List;

import ghidra.app.script.GhidraScript;
import ghidra.app.services.DataTypeManagerService;
import ghidra.framework.plugintool.PluginTool;
import ghidra.program.model.address.Address;
import ghidra.program.model.data.DataType;
import ghidra.program.model.data.DataTypeManager;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.Program;
import ghidra.program.model.symbol.SourceType;
import ghidra.program.model.symbol.Symbol;
import ghidra.program.model.symbol.SymbolTable;

public class ImportVarNames extends GhidraScript {
	
	private Program program;
	
	private DataType dt_word = null;
	private DataType dt_halfword = null;
	private DataType dt_byte = null;
	private DataType dt_void_p = null;
	
	private DataTypeManager getDataTypeManagerByName(String name) {
		//From FindDataTypeScript.java
		PluginTool tool = state.getTool();
		DataTypeManagerService service = tool.getService(DataTypeManagerService.class);
		DataTypeManager[] dataTypeManagers = service.getDataTypeManagers();
		for (DataTypeManager manager : dataTypeManagers) {
			String managerName = manager.getName();
			if (name.equals(managerName)) {
				return manager;
			}
		}
		return null;
	}
	
	private void findDataTypes() {
		DataTypeManager manager = getDataTypeManagerByName("BuiltInTypes");
		List<DataType> types = new LinkedList<DataType>();
		manager.getAllDataTypes(types);
		this.println("DataTypes: " + types.size());
		for(DataType dt : types) {
			String checkstr = dt.getName();
			//this.println("DataType checkstr: " + checkstr);
			if(checkstr.equalsIgnoreCase("dword")){ dt_word = dt;}
			else if(checkstr.equalsIgnoreCase("word")){ dt_halfword = dt;}
			else if(checkstr.equalsIgnoreCase("byte")) dt_byte = dt;
			else if(checkstr.equalsIgnoreCase("pointer")) dt_void_p = dt;
		}
	}
	
	protected void run() throws Exception {
		program = this.getCurrentProgram();
		if(program == null) {
			printerr("No active program...");
			return;
		}
		
		File infile = this.askFile("Select Input tsv", "Okay");
		//Expects NAME\tADDRESS\tTYPE\tARRAY_LEN
		
		findDataTypes();
		SymbolTable st = program.getSymbolTable();
		
		int count = 0;
		BufferedReader br = new BufferedReader(new FileReader(infile));
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.isBlank()) continue;
			if(line.startsWith("#")) continue;
			String[] split = line.split("\t");
			if(split.length < 2) continue;
			
			String sname = split[0];
			String araw = split[1];
			String eraw = split[3];
			int elementCount = 0;
			long fullSize = 0;
			DataType dt = dt_word;
			
			
			try {
				if(eraw.startsWith("0x")) {
					eraw = eraw.substring(2);
					elementCount = Integer.parseUnsignedInt(eraw.substring(2), 16);
				}
				else {
					elementCount = Integer.parseUnsignedInt(eraw);
				}
			}
			catch(NumberFormatException ex) {
				ex.printStackTrace();
				continue;
			}
			
			fullSize = elementCount << 2;
			if(split[2].endsWith("*")) {
				dt = dt_void_p;
			}
			else {
				//TODO char? ASCII strings?
				if(split[2].equalsIgnoreCase("u8") || split[2].equalsIgnoreCase("s8")) {
					dt = dt_byte;
					fullSize = elementCount;
				}
				if(split[2].equalsIgnoreCase("u16") || split[2].equalsIgnoreCase("s16")) {
					dt = dt_halfword;
					fullSize = elementCount << 1;
				}
			}
			
			Address myaddr = program.getAddressFactory().getAddress(araw);
			Address endaddr = myaddr.add(fullSize);

			//Clear existing symbols in region
			//clearListing(myaddr, endaddr);
			Symbol sym = null;
			Symbol[] symlist = st.getSymbols(myaddr);
			if((symlist != null) && (symlist.length > 0) && (symlist[0] != null)) {
				sym = symlist[0];
				//Change name
				sym.setName(sname, SourceType.DEFAULT);
			}
			else {
				sym = st.createLabel(myaddr, sname, SourceType.DEFAULT);
			}
			
			//Update data
			//TODO
			//How do I make an array???
			createData(myaddr, dt);
			
			count++;
		}
		br.close();
		
		println(count + " records processed!");
		
		
	}

}
